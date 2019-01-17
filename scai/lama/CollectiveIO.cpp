/**
 * @file CollectioveIO.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Parallel I/O operations for LAMA matrices and vectors using CollectiveFile.
 * @author Thomas Brandes
 * @date 20.06.2016
 */


#include <scai/lama/CollectiveIO.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/dmemo/CollectiveFile.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace lama
{

using dmemo::CollectiveFile;

using utilskernel::HArrayUtils;

using hmemo::_HArray;
using hmemo::HArray;

using common::ScalarType;

/** Internal id for DenseVector to verify correct file entries */

#define DENSE_VECTOR_CLASSID  1319816
#define CSR_MATRIX_CLASSID    1319827


/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CollectiveIO::logger, "CollectiveIO" )

/* --------------------------------------------------------------------------------- */

int CollectiveIO::getDenseVectorId()
{
    return DENSE_VECTOR_CLASSID;
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::write( CollectiveFile& file, const DenseVector<ValueType>& vector )
{
    const int val[] = { DENSE_VECTOR_CLASSID,
                        static_cast<int>( common::TypeTraits<IndexType>::stype ),
                        static_cast<int>( common::TypeTraits<ValueType>::stype )
                      };

    const auto& dist = vector.getDistribution();

    if ( dist.isReplicated() )
    {
        file.writeSingle( val, 3 );
        file.writeSingle( vector.size() );
        file.writeSingle( vector.getLocalValues() );
        return;
    }

    const auto& fileCommunicator   = file.getCommunicator();
    const auto& vectorCommunicator = dist.getCommunicator();

    SCAI_ASSERT_EQ_ERROR( fileCommunicator, vectorCommunicator, "serious communicator mismatch" )

    if ( dist.getBlockDistributionSize() == invalidIndex )
    {
        auto blockDist = dmemo::blockDistribution( vector.size(), file.getCommunicatorPtr() );
        DenseVector<ValueType> tmpVector;
        tmpVector.assignDistribute( vector, blockDist );
        write( file, tmpVector );
        return;
    }

    file.writeSingle( val, 3 );
    file.writeSingle( static_cast<IndexType>( vector.size() ) );
    file.writeAll( vector.getLocalValues() );
}

/* --------------------------------------------------------------------------------- */

template<typename TList>
struct IODataConverter;

template<>
struct IODataConverter<common::mepr::NullType>
{
    static void readAll( CollectiveFile&, _HArray&, const ScalarType stype, const IndexType, const IndexType )
    {
        COMMON_THROWEXCEPTION( "IODataConverter: " << stype << " unsupported" )
    }
};

template<typename ValueType, typename TailTypes>
struct IODataConverter<common::mepr::TypeList<ValueType, TailTypes> >
{
    static void readAll( CollectiveFile& file, _HArray& array, const ScalarType stype, const IndexType size, const IndexType offset )
    {
        if ( common::TypeTraits<ValueType>::stype == stype )
        {
            HArray<ValueType> typedArray;
            file.readAll( typedArray, size, offset );
            HArrayUtils::_assign( array, typedArray );
        }
        else
        {
            IODataConverter<TailTypes>::readAll( file, array, stype, size, offset );
        }
    }
};

/* --------------------------------------------------------------------------------- */

template<typename TList>
struct IOIndexConverter;

template<>
struct IOIndexConverter<common::mepr::NullType>
{
    static void readSingle( CollectiveFile&, IndexType&, const ScalarType stype )
    {
        COMMON_THROWEXCEPTION( "IOIndexConverter: " << stype << " unsupported" )
    }
};

template<typename ValueType, typename TailTypes>
struct IOIndexConverter<common::mepr::TypeList<ValueType, TailTypes> >
{
    static void readSingle( CollectiveFile& file, IndexType& value, const ScalarType stype )
    {
        if ( common::TypeTraits<ValueType>::stype == stype )
        {
            ValueType otherValue;
            file.readSingle( otherValue );
            value = static_cast<IndexType>( otherValue );
        }
        else
        {
            IOIndexConverter<TailTypes>::readSingle( file, value, stype );
        }
    }
};

/* --------------------------------------------------------------------------------- */

#define SCAI_INDEX_TYPE_LIST SCAI_TYPELIST( int, unsigned int, long, unsigned long )

template<typename ValueType>
void CollectiveIO::read( CollectiveFile& file, DenseVector<ValueType>& vector )
{
    int val[3];     // array to read the first three entries from the file

    file.readSingle( val, 3 );

    SCAI_ASSERT_EQ_ERROR( val[0], DENSE_VECTOR_CLASSID, "no dense vector in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType N;

    if ( fileIndexType == common::TypeTraits<IndexType>::stype )
    {
        file.readSingle( N );
    }
    else
    {
        IOIndexConverter<SCAI_INDEX_TYPE_LIST>::readSingle( file, N, fileIndexType );
    }

    auto dist = dmemo::blockDistribution( N, file.getCommunicatorPtr() );

    hmemo::HArray<ValueType> localValues;

    if ( fileDataType == common::TypeTraits<ValueType>::stype )
    {
        file.readAll( localValues, dist->getLocalSize(), dist->lb() );
    }
    else
    {
        IODataConverter<SCAI_ARRAY_TYPES_HOST_LIST>::readAll( file, localValues, fileDataType, dist->getLocalSize(), dist->lb() );
    }

    vector = DenseVector<ValueType>( dist, std::move( localValues ) );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeCSRMatrix( dmemo::CollectiveFile& file, const CSRSparseMatrix<ValueType>& csrMatrix )
{
    SCAI_ASSERT_ERROR( csrMatrix.getColDistribution().isReplicated(), "columns must be replicated." )

    SCAI_ASSERT_NE_ERROR( csrMatrix.getRowDistribution().getBlockDistributionSize(), invalidIndex, "no block distributon of rows" )

    const CSRStorage<ValueType>& csrLocal = csrMatrix.getLocalStorage();

    HArray<IndexType> csrSizes = csrLocal.getIA();

    HArrayUtils::unscan( csrSizes );  // we write sizes into file, not offsets

    const int val[] = { CSR_MATRIX_CLASSID,
                        static_cast<int>( common::TypeTraits<IndexType>::stype ),
                        static_cast<int>( common::TypeTraits<ValueType>::stype )
                      };

    file.writeSingle( val, 3 );
    file.writeSingle( csrMatrix.getNumRows() );
    file.writeSingle( csrMatrix.getNumColumns() );

    file.writeAll( csrSizes );
    file.writeAll( csrLocal.getJA() );
    file.writeAll( csrLocal.getValues() );
}
    
template<typename ValueType>
void CollectiveIO::write( dmemo::CollectiveFile& file, const SparseMatrix<ValueType>& matrix )
{
    if ( matrix.getFormat() == Format::CSR )
    {
        writeCSRMatrix( file, static_cast<const CSRSparseMatrix<ValueType>&>( matrix ) );
    }
    else
    {
        auto csrMatrix = convert<CSRSparseMatrix<ValueType>>( matrix );
        writeCSRMatrix( file, csrMatrix );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::readCSRMatrix( dmemo::CollectiveFile& file, CSRSparseMatrix<ValueType>& matrix )
{
    int val[3];     // array to read the first three entries from the file

    file.readSingle( val, 3 );

    SCAI_ASSERT_EQ_ERROR( val[0], CSR_MATRIX_CLASSID, "no CSR matrix in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType N, M;

    SCAI_ASSERT_EQ_ERROR( fileIndexType, common::TypeTraits<IndexType>::stype, "unsupported index type in input file." )
    SCAI_ASSERT_EQ_ERROR( fileDataType, common::TypeTraits<ValueType>::stype, "wrong value type in input file." )

    file.readSingle( N );
    file.readSingle( M );

    SCAI_LOG_ERROR( logger, file.getCommunicator() << ": read CSR matrix " << N << " x " << M )

    auto dist = dmemo::blockDistribution( N, file.getCommunicatorPtr() );
    auto localN = dist->getLocalSize();

    HArray<IndexType> csrIA( localN + 1 );

    file.readAll( csrIA, localN, dist->lb() );

    SCAI_LOG_ERROR( logger, file.getCommunicator() << ": read CSR IA " << csrIA )

    const IndexType localNNZ = HArrayUtils::scan1( csrIA );

    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    file.readAll( csrJA, localNNZ );
    SCAI_LOG_ERROR( logger, file.getCommunicator() << ": read CSR JA " << csrJA )
    file.readAll( csrValues, localNNZ );
    SCAI_LOG_ERROR( logger, file.getCommunicator() << ": read CSR Values " << csrValues )

    CSRStorage<ValueType> csrLocal( localN, M, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );
    matrix = CSRSparseMatrix<ValueType>( dist, std::move( csrLocal ) );
}

template<typename ValueType>
void CollectiveIO::read( dmemo::CollectiveFile& file, SparseMatrix<ValueType>& matrix )
{
    if ( matrix.getFormat() == Format::CSR )
    {
        readCSRMatrix( file, static_cast<CSRSparseMatrix<ValueType>&>( matrix ) );
    }
    else
    {
        CSRSparseMatrix<ValueType> tmpMatrix;
        readCSRMatrix( file, static_cast<CSRSparseMatrix<ValueType>&>( tmpMatrix ) );
        matrix = tmpMatrix;
    }
}

/* --------------------------------------------------------------------------------- */

#define SCAI_VECTOR_METHOD_INSTANTIATIONS( _type )          \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::write(                               \
        CollectiveFile& file,                               \
        const DenseVector<_type>& vector );                 \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::read(                                \
        CollectiveFile& file,                               \
        DenseVector<_type>& vector );                       \
                                                            \

SCAI_COMMON_LOOP( SCAI_VECTOR_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_VECTOR_METHOD_INSTANTIATIONS

#define SCAI_MATRIX_METHOD_INSTANTIATIONS( _type )          \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::write(                               \
        CollectiveFile& file,                               \
        const SparseMatrix<_type>& vector );                \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::read(                                \
        CollectiveFile& file,                               \
        SparseMatrix<_type>& vector );                      \

SCAI_COMMON_LOOP( SCAI_MATRIX_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_MATRIX_METHOD_INSTANTIATIONS

}  // namespace lama

}  // namespace scai
