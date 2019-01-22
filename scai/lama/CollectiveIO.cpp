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
#include <scai/lama/SparseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/dmemo/CollectiveFile.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>

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

#define DENSE_VECTOR_CLASSID  0x4711E00
#define SPARSE_VECTOR_CLASSID 0x4711E01
#define DENSE_MATRIX_CLASSID  0x4711E10
#define CSR_MATRIX_CLASSID    0x4711E11


/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( CollectiveIO::logger, "CollectiveIO" )

/* --------------------------------------------------------------------------------- */

int CollectiveIO::getDenseVectorId()
{
    return DENSE_VECTOR_CLASSID;
}

/* --------------------------------------------------------------------------------- */
/*   write Vector                                                                    */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::write( 
    CollectiveFile& file, 
    const Vector<ValueType>& vector,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType )
{
    auto kind = vector.getVectorKind();

    if ( kind == VectorKind::DENSE )
    {
        writeDenseVector( file, static_cast<const DenseVector<ValueType>&>( vector ), fileIndexType, fileDataType );
    }
    else if ( kind == VectorKind::SPARSE )
    {
        writeSparseVector( file, static_cast<const SparseVector<ValueType>&>( vector ), fileIndexType, fileDataType );
    }
    else 
    {
        COMMON_THROWEXCEPTION( "Unsupported vector kind: " << kind );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::write( 
    CollectiveFile& file, 
    const Vector<ValueType>& vector )
{
    write( file, vector, common::TypeTraits<IndexType>::stype, common::TypeTraits<ValueType>::stype );
}

/* --------------------------------------------------------------------------------- */
/*   write Matrix                                                                    */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::write( 
    CollectiveFile& file, 
    const Matrix<ValueType>& matrix,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType )
{
    auto kind = matrix.getMatrixKind();

    if ( kind == MatrixKind::DENSE )
    {
        writeDenseMatrix( file, static_cast<const DenseMatrix<ValueType>&>( matrix ), fileIndexType, fileDataType );
    }
    else if ( kind == MatrixKind::SPARSE )
    {
        writeSparseMatrix( file, static_cast<const SparseMatrix<ValueType>&>( matrix ), fileIndexType, fileDataType );
    }
    else 
    {
        COMMON_THROWEXCEPTION( "Unsupported matrix kind: " << kind );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::write( 
    CollectiveFile& file, 
    const Matrix<ValueType>& matrix )
{
    write( file, matrix, common::TypeTraits<IndexType>::stype, common::TypeTraits<ValueType>::stype );
}

/* --------------------------------------------------------------------------------- */
    
template<typename ValueType>
void CollectiveIO::writeDenseVector( 
    CollectiveFile& file, 
    const DenseVector<ValueType>& vector,
    const common::ScalarType fileIndexType, 
    const common::ScalarType fileDataType )
{
    const auto& dist = vector.getDistribution();

    const auto& fileCommunicator   = file.getCommunicator();
    const auto& vectorCommunicator = dist.getCommunicator();

    SCAI_ASSERT_EQ_ERROR( fileCommunicator, vectorCommunicator, "serious communicator mismatch" )

    if ( dist.getBlockDistributionSize() == invalidIndex || fileCommunicator != vectorCommunicator )
    {
        auto blockDist = dmemo::blockDistribution( vector.size(), file.getCommunicatorPtr() );
        DenseVector<ValueType> tmpVector;
        tmpVector.assignDistribute( vector, blockDist );
        writeIt( file, tmpVector, fileIndexType, fileDataType );
    }
    else
    {
        writeIt( file, vector, fileIndexType, fileDataType );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeSparseVector(
    CollectiveFile& file,
    const SparseVector<ValueType>& vector,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType )
{
    const auto& dist = vector.getDistribution();

    const auto& fileCommunicator   = file.getCommunicator();
    const auto& vectorCommunicator = dist.getCommunicator();

    if ( dist.getBlockDistributionSize() == invalidIndex || fileCommunicator != vectorCommunicator )
    {
        auto blockDist = dmemo::blockDistribution( vector.size(), file.getCommunicatorPtr() );

        SCAI_LOG_WARN( logger, "redistribute sparse vector : " << vector << " for collective IO with new dist " << *blockDist )

        SparseVector<ValueType> tmpVector;
        tmpVector.assignDistribute( vector, blockDist );
        writeIt( file, tmpVector, fileIndexType, fileDataType );
    }
    else
    {
        writeIt( file, vector, fileIndexType, fileDataType );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeDenseMatrix( 
    CollectiveFile& file, 
    const DenseMatrix<ValueType>& matrix,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType  )
{
    // writing into file only possible for block distributed matrices 

    const auto& rowDist = matrix.getRowDistribution();
    const auto& colDist = matrix.getColDistribution();

    if ( rowDist.getBlockDistributionSize() == invalidIndex || file.getCommunicator() != rowDist.getCommunicator() )
    {
        // block distribution onto processor set of file communicator required

        auto blockDist = dmemo::blockDistribution( matrix.getNumRows(), file.getCommunicatorPtr() );
        auto repDist = dmemo::noDistribution( matrix.getNumColumns() );
        auto tmpMatrix = distribute<DenseMatrix<ValueType>>( matrix, blockDist, repDist );
        writeIt( file, tmpMatrix, fileIndexType, fileDataType );
    } 
    else if ( !colDist.isReplicated() )
    {
        auto blockDist = matrix.getRowDistributionPtr();
        auto repDist = dmemo::noDistribution( matrix.getNumColumns() );
        auto tmpMatrix = distribute<DenseMatrix<ValueType>>( matrix, blockDist, repDist );
        writeIt( file, tmpMatrix, fileIndexType, fileDataType );
    } 
    else
    {
        writeIt( file, matrix, fileIndexType, fileDataType );
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeSparseMatrix(
    CollectiveFile& file,
    const SparseMatrix<ValueType>& matrix,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType  )
{
    // writing into file only possible for block distributed matrices 

    const auto& rowDist = matrix.getRowDistribution();
    const auto& colDist = matrix.getColDistribution();

    if ( rowDist.getBlockDistributionSize() == invalidIndex || file.getCommunicator() != rowDist.getCommunicator() )
    {
        // block distribution onto processor set of file communicator required

        auto blockDist = dmemo::blockDistribution( matrix.getNumRows(), file.getCommunicatorPtr() );
        auto repDist = dmemo::noDistribution( matrix.getNumColumns() );
        auto tmpMatrix = distribute<CSRSparseMatrix<ValueType>>( matrix, blockDist, repDist );
        writeIt( file, tmpMatrix, fileIndexType, fileDataType );
    }
    else if ( !colDist.isReplicated() )
    {
        auto blockDist = matrix.getRowDistributionPtr();
        auto repDist = dmemo::noDistribution( matrix.getNumColumns() );
        auto tmpMatrix = distribute<CSRSparseMatrix<ValueType>>( matrix, blockDist, repDist );
        writeIt( file, tmpMatrix, fileIndexType, fileDataType );
    }
    else if ( matrix.getFormat() != Format::CSR )
    {
        auto tmpMatrix = convert<CSRSparseMatrix<ValueType>>( matrix );
        writeIt( file, tmpMatrix, fileIndexType, fileDataType );
    }
    else 
    {
        writeIt( file, static_cast<const CSRSparseMatrix<ValueType>&>( matrix ), fileIndexType, fileDataType );
    }
}

/* --------------------------------------------------------------------------------- */
/*   write/read dense vector                                                         */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeIt( 
    CollectiveFile& file, 
    const DenseVector<ValueType>& denseVector,
    const common::ScalarType fileIndexType, 
    const common::ScalarType fileDataType )
{
    const auto& dist = denseVector.getDistribution();

    SCAI_ASSERT_EQ_ERROR( file.getCommunicator(), dist.getCommunicator(), "serious mismatch" )
    SCAI_ASSERT_NE_ERROR( dist.getBlockDistributionSize(), invalidIndex, "no block distribution" )

    const int val[] = { DENSE_VECTOR_CLASSID,
                        static_cast<int>( fileIndexType ),
                        static_cast<int>( fileDataType )
                      };

    file.writeSingle( val, 3 );
    file.writeSingle( denseVector.size(), fileIndexType );
    file.writeAll( denseVector.getLocalValues(), fileDataType );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::readIt( CollectiveFile& file, DenseVector<ValueType>& denseVector )
{
    int val[3];     // array to read the first three entries from the file

    file.readSingle( val, 3 );

    SCAI_ASSERT_EQ_ERROR( val[0], DENSE_VECTOR_CLASSID, "no dense vector in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType N;

    file.readSingle( N, fileIndexType );

    auto dist = dmemo::blockDistribution( N, file.getCommunicatorPtr() );

    hmemo::HArray<ValueType> localValues;

    file.readAll( localValues, dist->getLocalSize(), dist->lb(), fileDataType );

    denseVector = DenseVector<ValueType>( dist, std::move( localValues ) );
}

/* --------------------------------------------------------------------------------- */
/*   write/read sparse vector                                                        */
/*   <SPARSE_VECTOR_CLASSID><indexType><dataType><N><NNZ><zero>                      */
/*   <[non-zero-indexes]><[non-zero-values]>                                         */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeIt(
    CollectiveFile& file,
    const SparseVector<ValueType>& sparseVector,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType )
{
    const auto& dist = sparseVector.getDistribution();

    SCAI_ASSERT_EQ_ERROR( file.getCommunicator(), dist.getCommunicator(), "serious mismatch" )
    SCAI_ASSERT_NE_ERROR( dist.getBlockDistributionSize(), invalidIndex, "no block distribution" )

    const int val[] = { SPARSE_VECTOR_CLASSID,
                        static_cast<int>( fileIndexType ),
                        static_cast<int>( fileDataType )
                      };

    file.writeSingle( val, 3 );

    HArray<IndexType> nonZeroIndexes;   // global non-zero indexes
    dist.local2GlobalV( nonZeroIndexes, sparseVector.getNonZeroIndexes() );
    const HArray<ValueType>& nonZeroValues = sparseVector.getNonZeroValues();

    SCAI_ASSERT_EQ_ERROR( nonZeroIndexes.size(), nonZeroValues.size(), "Inconsistent sparse vector: " << sparseVector );

    const IndexType N   = sparseVector.size();
    const IndexType NNZ = file.getCommunicator().sum( nonZeroIndexes.size() );

    file.writeSingle( N, fileIndexType );
    file.writeSingle( NNZ, fileIndexType );

    file.writeSingle( sparseVector.getZero(), fileDataType );

    SCAI_LOG_TRACE( logger, file.getCommunicator() << ": write indexes = " << printIt( nonZeroIndexes )
                             << ": write values = " << printIt( nonZeroValues ) )

    file.writeAll( nonZeroIndexes, fileIndexType );
    file.writeAll( nonZeroValues, fileDataType );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::readIt( CollectiveFile& file, SparseVector<ValueType>& sparseVector )
{
    int val[3];     // array to read the first three entries from the file

    file.readSingle( val, 3 );

    SCAI_ASSERT_EQ_ERROR( val[0], SPARSE_VECTOR_CLASSID, "no sparse vector in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType N;     // size of the vector
    IndexType NNZ;   // number of non-zeros in the vector
    ValueType zero;  // zero value

    file.readSingle( N, fileIndexType );
    file.readSingle( NNZ, fileIndexType );
    file.readSingle( zero, fileDataType );

    SCAI_LOG_INFO( logger, file.getCommunicator() << ": read sparse vector: N = " << N << ", NNZ = " << NNZ << ", zero = " << zero )

    // parallel read of the non-zero indexes + values

    auto dist = dmemo::blockDistribution( NNZ, file.getCommunicatorPtr() );

    hmemo::HArray<IndexType> nonZeroIndexes;
    hmemo::HArray<ValueType> nonZeroValues;

    file.readAll( nonZeroIndexes, dist->getLocalSize(), dist->lb(), fileIndexType );
    file.readAll( nonZeroValues, dist->getLocalSize(), dist->lb(), fileDataType );

    SCAI_LOG_TRACE( logger, file.getCommunicator() << ": read indexes = " << hmemo::printIt( nonZeroIndexes ) 
                            << ": read values = " << hmemo::printIt( nonZeroValues ) )

    IndexType offset = nonZeroIndexes.size() > 0 ? nonZeroIndexes[0] : invalidIndex;

    auto newDist = dmemo::genBlockDistributionByOffset( N, offset, file.getCommunicatorPtr() );

    newDist->global2LocalV( nonZeroIndexes, nonZeroIndexes );

    sparseVector = SparseVector<ValueType>( newDist, std::move( nonZeroIndexes ), std::move( nonZeroValues ), zero );
}

/* --------------------------------------------------------------------------------- */
/*   write/read CSR sparse matrix                                                    */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeIt( 
    dmemo::CollectiveFile& file, 
    const CSRSparseMatrix<ValueType>& csrMatrix,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType )
{
    const CSRStorage<ValueType>& csrLocal = csrMatrix.getLocalStorage();

    HArray<IndexType> csrSizes = csrLocal.getIA();

    HArrayUtils::unscan( csrSizes );  // we write sizes into file, not offsets

    const int val[] = { CSR_MATRIX_CLASSID,
                        static_cast<int>( fileIndexType ),
                        static_cast<int>( fileDataType )
                      };

    file.writeSingle( val, 3 );
    file.writeSingle( csrMatrix.getNumRows(), fileIndexType );
    file.writeSingle( csrMatrix.getNumColumns(), fileIndexType );

    file.writeAll( csrSizes, fileIndexType );
    file.writeAll( csrLocal.getJA(), fileIndexType );
    file.writeAll( csrLocal.getValues(), fileDataType );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::readIt( dmemo::CollectiveFile& file, CSRSparseMatrix<ValueType>& matrix )
{
    int val[3];     // array to read the first three entries from the file

    file.readSingle( val, 3 );

    SCAI_LOG_INFO( logger, "read header for CSR matrix: " << val[0] << ", " << val[1] << ", " << val[2] )

    SCAI_ASSERT_EQ_ERROR( val[0], CSR_MATRIX_CLASSID, "no CSR matrix in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType N, M;

    file.readSingle( N, fileIndexType );
    file.readSingle( M, fileIndexType );

    SCAI_LOG_DEBUG( logger, file.getCommunicator() << ": read CSR matrix " << N << " x " << M )

    auto dist = dmemo::blockDistribution( N, file.getCommunicatorPtr() );
    auto localN = dist->getLocalSize();

    HArray<IndexType> csrIA( localN + 1 );

    file.readAll( csrIA, localN, dist->lb(), fileIndexType );

    const IndexType localNNZ = HArrayUtils::scan1( csrIA );

    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    file.readAll( csrJA, localNNZ, fileIndexType );
    file.readAll( csrValues, localNNZ, fileDataType );

    CSRStorage<ValueType> csrLocal( localN, M, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );
    matrix = CSRSparseMatrix<ValueType>( dist, std::move( csrLocal ) );
}

/* --------------------------------------------------------------------------------- */
/*   write/read dense matrix                                                         */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::writeIt( 
    CollectiveFile& file, 
    const DenseMatrix<ValueType>& matrix,
    const common::ScalarType fileIndexType,
    const common::ScalarType fileDataType  )
{
    // matrix is block distributed onto the processors of the file communicator

    const int val[] = { DENSE_MATRIX_CLASSID,
                        static_cast<int>( fileIndexType ),
                        static_cast<int>( fileDataType )
                      };

    file.writeSingle( val, 3 );
   
    file.writeSingle( matrix.getNumRows(), fileIndexType );
    file.writeSingle( matrix.getNumColumns(), fileIndexType );

    const HArray<ValueType>& denseData = matrix.getLocalStorage().getValues();

    file.writeAll( denseData, fileDataType );
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::readIt( dmemo::CollectiveFile& file, DenseMatrix<ValueType>& matrix )
{
    int val[3];     // array to read the first three entries from the file

    file.readSingle( val, 3 );

    SCAI_ASSERT_EQ_ERROR( val[0], DENSE_MATRIX_CLASSID, "no DENSE matrix in file" )

    auto fileIndexType = ScalarType( val[1] );
    auto fileDataType  = ScalarType( val[2] );

    IndexType N, M;

    file.readSingle( N, fileIndexType );
    file.readSingle( M, fileIndexType );

    SCAI_LOG_INFO( logger, file.getCommunicator() << ": read CSR matrix " << N << " x " << M )

    auto dist = dmemo::blockDistribution( N, file.getCommunicatorPtr() );
    auto localN = dist->getLocalSize();

    HArray<ValueType> denseValues;

    file.readAll( denseValues, localN * M, fileDataType );

    DenseStorage<ValueType> denseLocal( localN, M, std::move( denseValues ) );
    matrix = DenseMatrix<ValueType>( dist, std::move( denseLocal ) );
}

/* --------------------------------------------------------------------------------- */
/*  read  Vector / Matrix                                                            */
/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::read( CollectiveFile& file, Vector<ValueType>& vector )
{
    auto kind = vector.getVectorKind();

    if ( kind == VectorKind::SPARSE ) 
    {
        readIt( file, static_cast<SparseVector<ValueType>&>( vector ) );
    }
    else if ( kind == VectorKind::DENSE ) 
    {
        readIt( file, static_cast<DenseVector<ValueType>&>( vector ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "unsupported vector kind " << kind << " for collective read" )
    }
}

/* --------------------------------------------------------------------------------- */

template<typename ValueType>
void CollectiveIO::read( dmemo::CollectiveFile& file, Matrix<ValueType>& matrix )
{
    if ( matrix.getFormat() == Format::CSR )
    {
        readIt( file, static_cast<CSRSparseMatrix<ValueType>&>( matrix ) );
    }
    else if ( matrix.getFormat() == Format::DENSE )
    {
        readIt( file, static_cast<DenseMatrix<ValueType>&>( matrix ) );
    }
    else
    {
        CSRSparseMatrix<ValueType> tmpMatrix;
        readIt( file, static_cast<CSRSparseMatrix<ValueType>&>( tmpMatrix ) );
        matrix = tmpMatrix;
    }
}

/* --------------------------------------------------------------------------------- */
/*  Method instantiations                                                            */
/* --------------------------------------------------------------------------------- */

#define SCAI_VECTOR_METHOD_INSTANTIATIONS( _type )          \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::write(                               \
        CollectiveFile& file,                               \
        const Vector<_type>& vector );                      \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::write(                               \
        CollectiveFile& file,                               \
        const Vector<_type>& vector,                        \
        const common::ScalarType,                           \
        const common::ScalarType );                         \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::read(                                \
        CollectiveFile& file,                               \
        Vector<_type>& vector );                            \
                                                            \

SCAI_COMMON_LOOP( SCAI_VECTOR_METHOD_INSTANTIATIONS, SCAI_ARRAY_TYPES_HOST )

#undef SCAI_VECTOR_METHOD_INSTANTIATIONS

#define SCAI_MATRIX_METHOD_INSTANTIATIONS( _type )          \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::write(                               \
        CollectiveFile& file,                               \
        const Matrix<_type>& vector );                      \
                                                            \
    template COMMON_DLL_IMPORTEXPORT                        \
    void CollectiveIO::read(                                \
        CollectiveFile& file,                               \
        Matrix<_type>& vector );                            \

SCAI_COMMON_LOOP( SCAI_MATRIX_METHOD_INSTANTIATIONS, SCAI_NUMERIC_TYPES_HOST )

#undef SCAI_MATRIX_METHOD_INSTANTIATIONS

}  // namespace lama

}  // namespace scai
