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

#include <scai/dmemo/CollectiveFile.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/common/TypeTraits.hpp>

namespace scai
{

namespace lama
{

using dmemo::CollectiveFile;

using hmemo::_HArray;
using hmemo::HArray;

using common::ScalarType;

/** Internal id for DenseVector to verify correct file entries */

#define DENSE_VECTOR_CLASSID 1319816

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
            utilskernel::HArrayUtils::_assign( array, typedArray );
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

#define SCAI_COL_IO_METHOD_INSTANTIATIONS( _type )          \
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


SCAI_COMMON_LOOP( SCAI_COL_IO_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

}  // namespace lama

}  // namespace scai
