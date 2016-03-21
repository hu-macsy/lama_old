/**
 * @file SparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Template specilization of the matrix template for distributed matrixes.
 * @author Jiri Kraus, Thomas Brandes
 * @date 02.04.2012
 * @since 1.0.0
 */

// hpp
#include <scai/lama/matrix/SparseMatrix.hpp>

// local library

#include <scai/lama/matrix/DenseMatrix.hpp>

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>


// scai internal libraries
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
 
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/macros/throw.hpp>
#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <cmath>

using namespace scai::hmemo;
using namespace scai::dmemo;

namespace scai
{

namespace lama
{

using common::shared_ptr;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using sparsekernel::CSRKernelTrait;

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SparseMatrix<ValueType>::logger, "Matrix.SparseMatrix" )

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( common::shared_ptr<MatrixStorage<ValueType> > storage )
    :

    CRTPMatrix<SparseMatrix<ValueType>,ValueType>( storage->getNumRows(), storage->getNumColumns() )
{
    mLocalData = storage;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( storage->newMatrixStorage() );
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( common::shared_ptr<MatrixStorage<ValueType> > storage, DistributionPtr rowDist )
    :

    CRTPMatrix<SparseMatrix<ValueType>,ValueType>(
       rowDist, DistributionPtr( new NoDistribution( storage->getNumColumns() ) ) )
{
    mLocalData = storage;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( storage->newMatrixStorage() );

    if( storage->getNumRows() == rowDist->getLocalSize() )
    {
        // data fits, no more to do
    }
    else if( storage->getNumRows() == rowDist->getGlobalSize() )
    {
        // localize the data according to row distribution, use splitHalo with replicated columns

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), rowDist.get() );
    }
    else
    {
        COMMON_THROWEXCEPTION( *storage << ": does not fit to row distribution " << *rowDist )
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix(
    common::shared_ptr<MatrixStorage<ValueType> > localData,
    DistributionPtr rowDist,
    DistributionPtr colDist )
    :

    CRTPMatrix<SparseMatrix<ValueType>,ValueType>( rowDist, colDist )
{
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getGlobalSize() )

    mLocalData = localData;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( localData->newMatrixStorage() );

    if( localData->getNumRows() == rowDist->getLocalSize() )
    {
        // build just the halo

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );

    }
    else if( localData->getNumRows() == rowDist->getGlobalSize() )
    {
        // we also localize the data

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, rowDist.get() );
    }
    else
    {
        COMMON_THROWEXCEPTION( *localData << ": does not fit to row distribution " << *rowDist )
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix(
    common::shared_ptr<MatrixStorage<ValueType> > localData,
    common::shared_ptr<MatrixStorage<ValueType> > haloData,
    const Halo& halo,
    DistributionPtr rowDist,
    DistributionPtr colDist )
    :

    CRTPMatrix<SparseMatrix<ValueType>,ValueType>( rowDist, colDist )
{
    SCAI_LOG_INFO( logger, "Construct sparse matrix with finalized local, halo storage + Halo" )

    // TODO: asserts for correct sizes of all relevant sizes

    SCAI_ASSERT_EQUAL_ERROR( localData->getNumRows(), rowDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getLocalSize() )

    SCAI_ASSERT_EQUAL_ERROR( haloData->getNumRows(), rowDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( haloData->getNumColumns(), halo.getHaloSize() )

    // done by constructor for CRTPMatrix<SparseMatrix<ValueType>, ValueType >:  Matrix::setDistributedMatrix( rowDist, colDist );

    mLocalData = localData; // flat copy
    mHaloData = haloData; // flat copy
    mHalo = halo;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=( const SparseMatrix& matrix )
{
    SCAI_LOG_INFO( logger, " = SparseMatrix : " << matrix )
    assign( matrix );
    return *this;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void SparseMatrix<ValueType>::checkSettings()
{
    Matrix::checkSettings();

    if( mHaloData->getNumRows() )
    {
        // If halo storage is available, its size must fit with local storage

        SCAI_ASSERT_EQUAL_DEBUG( mLocalData->getNumRows(), mHaloData->getNumRows() )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseMatrix<ValueType>::isConsistent() const
{
    int consistencyErrors = 0;

    // ToDo: this implementation should use a corresponding predicate of MatrixStorage

    try
    {
        Matrix::checkSettings();

        SCAI_ASSERT_EQUAL_ERROR( getRowDistribution().getLocalSize(), mLocalData->getNumRows() )
        SCAI_ASSERT_EQUAL_ERROR( mHaloData->getNumRows(), mLocalData->getNumRows() )

        mLocalData->check( "check for consistency" );
        mHaloData->check( "check for consistency" );

        // ToDo: check Halo
    }
    catch( ... )
    {
        consistencyErrors = 1;
    }

    // use communicator for global reduction to make sure that all processors return same value.

    consistencyErrors = getRowDistribution().getCommunicator().sum( consistencyErrors );

    return 0 == consistencyErrors;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( const Matrix& other, const bool transposeFlag /* = false */)

{
    Matrix::inheritAttributes( other ); // context, communication, compute kind

    if( transposeFlag )
    {
        assignTranspose( other );
    }
    else
    {
        assign( other );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( const Matrix& other, DistributionPtr rowDist, DistributionPtr colDist )

{
    Matrix::inheritAttributes( other ); // context, communication, compute kind

    assign( other );
    this->redistribute( rowDist, colDist );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( const SparseMatrix<ValueType>& other )
    :

    CRTPMatrix<SparseMatrix<ValueType>,ValueType>( other )

{
    // instead of calling assign( other ), we copy directly to avoid type queries

    // make deep copies of the storage data

    mLocalData = shared_ptr<MatrixStorage<ValueType> >( other.getLocalStorage().copy() );
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( other.getHaloStorage().copy() );

    // just copy the halo

    mHalo = other.getHalo();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;

    mLocalData->clear();
    mHaloData->clear();
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::purge()
{
    // Note: purge will free the memory and not only reset sizes

    mNumRows = 0;
    mNumColumns = 0;

    mLocalData->purge();
    mHaloData->purge();
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::allocate( DistributionPtr distribution, DistributionPtr colDistribution )
{
    SCAI_LOG_DEBUG( logger,
                    *this << ": allocate with row dist = " << *distribution << ", col dist = " << *colDistribution )

    Matrix::setDistributedMatrix( distribution, colDistribution );

    IndexType localSize = distribution->getLocalSize();
    IndexType localColSize = colDistribution->getLocalSize();

    mLocalData->allocate( localSize, localColSize );
    mHaloData->allocate( localSize, colDistribution->getGlobalSize() );
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::allocate( const IndexType numRows, const IndexType numColumns )
{
    Matrix::setReplicatedMatrix( numRows, numColumns );

    mLocalData->allocate( numRows, numColumns );
    mHaloData->allocate( numRows, 0 );
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const Matrix& matrix )
{
    SCAI_LOG_INFO( logger, "assign " << matrix << " to " << *this )

    this->setContextPtr( matrix.getContextPtr() );

    const SparseMatrix<ValueType>* sparseMatrix = dynamic_cast<const SparseMatrix<ValueType>*>( &matrix );

    if( sparseMatrix )
    {
        // for a sparse matrix local + halo part can be assigned

        assign( *sparseMatrix );
    }
    else
    {
        // convert dense matrix to a sparse matrix

        DistributionPtr colDist = matrix.getColDistributionPtr();
        DistributionPtr rowDist = matrix.getRowDistributionPtr();

        Matrix::setDistributedMatrix( rowDist, colDist );

        matrix.buildLocalStorage( *mLocalData ); // local storage with all columns

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assignTranspose( const Matrix& matrix )
{
    SCAI_LOG_INFO( logger, "assign transposed " << matrix << " to " << *this )

    const SparseMatrix<ValueType>* sparseMatrix = dynamic_cast<const SparseMatrix<ValueType>*>( &matrix );

    if( sparseMatrix )
    {
        assignTransposeImpl( *sparseMatrix );
    }
    else
    {
        COMMON_THROWEXCEPTION( "SparseMatrix::assign currently only implemented for sparse matrices of same type" )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assignTransposeImpl( const SparseMatrix<ValueType>& matrix )

{
    SCAI_LOG_INFO( logger, "transpose sparse matrix with same value type, switch row/col distributions" )

    // assign matrix properties

    Matrix::setDistributedMatrix( matrix.getColDistributionPtr(), matrix.getRowDistributionPtr() );

    if( getRowDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, "transpose local storage, input = " << matrix.getLocalStorage() )
        mLocalData->assignTranspose( matrix.getLocalStorage() );
        SCAI_LOG_DEBUG( logger, "transposed local storage, is = " << *mLocalData )
        mHaloData->allocate( getRowDistribution().getLocalSize(), 0 );
        mHalo.clear();
        SCAI_LOG_DEBUG( logger, "transposed halo storage, is = " << *mHaloData )
    }
    else if( getRowDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "transpose not supported for replicated matrices with distributed columns" )
    }
    else if( getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "transpose not supported for distributed matrices with replicated columns" )
    }
    else
    {
        // rows and columns are both distributed

        SCAI_LOG_INFO( logger, "local transpose of " << matrix.getLocalStorage() )

        mLocalData->assignTranspose( matrix.getLocalStorage() );

        SCAI_LOG_INFO( logger, "local transposed = " << mLocalData )

        SCAI_LOG_INFO( logger, "halo transpose of " << matrix.getHaloStorage() )

        HArray<IndexType> sendIA;
        HArray<IndexType> sendSizes;
        HArray<IndexType> sendJA;
        HArray<ValueType> sendValues;

        matrix.getHaloStorage().buildCSCData( sendIA, sendJA, sendValues );

        SCAI_LOG_DEBUG( logger,
                        matrix.getHaloStorage() << ": CSC data, IA = " << sendIA << ", JA = " << sendJA << ", Values = " << sendValues )

        // for initial communication we need the sizes and not the offsets

        _MatrixStorage::offsets2sizes( sendSizes, sendIA );

        SCAI_LOG_DEBUG( logger, "sendSizes with " << sendSizes.size() << " entries" )

        HArray<IndexType> recvSizes;

        const Halo& matrixHalo = matrix.getHalo();

        SCAI_ASSERT_EQUAL_DEBUG( sendSizes.size(), matrix.getHaloStorage().getNumColumns() )

        // send all other processors the number of columns

        const Communicator& comm = getRowDistribution().getCommunicator();

        // Use communication plans of halo in inverse manner to send col sizes

        const CommunicationPlan& sendSizesPlan = matrixHalo.getRequiredPlan();
        const CommunicationPlan& recvSizesPlan = matrixHalo.getProvidesPlan();

        comm.exchangeByPlan( recvSizes, recvSizesPlan, sendSizes, sendSizesPlan );

        ContextPtr contextPtr = Context::getHostPtr();

        // Now we know the sizes, we can pack the data

        HArray<IndexType> recvJA;
        HArray<ValueType> recvValues;

        // Before send of JA: translate back the local (row) indexes to global indexes

        {
            const Distribution& dist = matrix.getRowDistribution();
            const IndexType nJA = sendJA.size();

            WriteAccess<IndexType> ja( sendJA, contextPtr );

            for( IndexType jj = 0; jj < nJA; ++jj )
            {
                ja[jj] = dist.local2global( ja[jj] );
            }
        }

        // Build the communication plan for sparse data

        {
            ReadAccess<IndexType> sendColSizes( sendSizes, contextPtr );
            ReadAccess<IndexType> recvColSizes( recvSizes, contextPtr );

            CommunicationPlan sendDataPlan( sendSizesPlan, sendColSizes.get() );
            CommunicationPlan recvDataPlan( recvSizesPlan, recvColSizes.get() );

            comm.exchangeByPlan( recvJA, recvDataPlan, sendJA, sendDataPlan );
            comm.exchangeByPlan( recvValues, recvDataPlan, sendValues, sendDataPlan );
        }

        // Now we have to build the halo by merging the rows of all received data

        HArray<IndexType> haloRowSizes;
        HArray<IndexType> haloJA;
        HArray<ValueType> haloValues;

        const HArray<IndexType>& rowIndexes = matrixHalo.getProvidesIndexes();

        const IndexType mNumLocalRows = getRowDistribution().getLocalSize();

        {
            // allocate rowSizes as offset array to avoid reallocation
            WriteOnlyAccess<IndexType> rowSizes( haloRowSizes, mNumLocalRows + 1 );
        }

        MatrixStorage<ValueType>::joinRows( haloRowSizes, haloJA, haloValues, mNumLocalRows, rowIndexes, recvSizes,
                                            recvJA, recvValues );

        // Now this partition has all sparse data from other partitions to build my new halo

        mHaloData->setCompressThreshold( 1.0f ); // halo rows might be compressed

        _MatrixStorage::sizes2offsets( haloRowSizes );

        mHaloData->setCSRData( mNumLocalRows, mNumColumns, haloJA.size(), haloRowSizes, haloJA, haloValues );

        // Now build a new halo and localize columns in mHaloData

        mHaloData->buildHalo( mHalo, getColDistribution() );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const SparseMatrix<ValueType>& matrix )
{
    if( this == &matrix )
    {
        SCAI_LOG_INFO( logger, "self assign sparse matrix = " << matrix )
    }
    else
    {
        SCAI_LOG_INFO( logger, "copy/convert assign sparse matrix = " << matrix )
    }

    Matrix::setDistributedMatrix( matrix.getRowDistributionPtr(), matrix.getColDistributionPtr() );

    mLocalData->assign( matrix.getLocalStorage() );

    SCAI_LOG_DEBUG( logger, "assigned local storage, my local = " << *mLocalData )

    const MatrixStorage<ValueType>&  matrixHaloData = matrix.getHaloStorage();

    if (     ( matrixHaloData.getNumRows() > 0 )
          || ( matrixHaloData.getNumColumns()  > 0 ) )
    {
        mHaloData->assign( matrixHaloData );
    }
    else
    {
        // we set the halo size to 0 x 0 instead of n x 0 to avoid allocation of ia data

        mHaloData->clear();   
    }

    mHalo = matrix.getHalo();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const _MatrixStorage& storage )
{
    SCAI_LOG_INFO( logger, "assign matrix storage = " << storage )

    const IndexType numRows = storage.getNumRows();
    const IndexType numColumns = storage.getNumColumns();

    Matrix::setReplicatedMatrix( numRows, numColumns );

    // TODO: allow flexibility regarding the context, e.g. format conversion should be done on GPU

    mLocalData->assign( storage );
    mHaloData->allocate( numRows, 0 ); // empty halo storage, halo
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist )
{
    SCAI_LOG_INFO( logger,
                   "assign storage = " << storage << ", row dist = " << *rowDist << ", col dist = " << *colDist );

    SCAI_ASSERT_EQUAL_ERROR( storage.getNumColumns(), colDist->getGlobalSize() )

    // split can only be applied to storage of same type, so check if we need to convert value type

    const MatrixStorage<ValueType>* typedStorage = dynamic_cast<const MatrixStorage<ValueType>*>( &storage );

    if( !typedStorage )
    {
        mLocalData->assign( storage ); // conversion
        typedStorage = mLocalData.get();
    }

    Matrix::setDistributedMatrix( rowDist, colDist );

    const Communicator& comm = rowDist->getCommunicator();

    bool isLocal = storage.getNumRows() == rowDist->getLocalSize();
    bool isGlobal = storage.getNumRows() == rowDist->getGlobalSize();

    bool allLocal = comm.all( isLocal ); // local storage on all processors
    bool allGlobal = comm.all( isGlobal ); // global storage on all processors

    if( allLocal )
    {
        // each processor has exactly its local part

        typedStorage->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }

    else if( allGlobal )
    {
        // we also localize the rows of the matrix

        typedStorage->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, rowDist.get() );
    }

    else
    {
        // Note: all processors with throw this exception

        COMMON_THROWEXCEPTION(
            comm << ": " << storage << ": does not fit to row distribution " << *rowDist << ", isLocal = " << isLocal << ", isGlobal = " << isGlobal )
    }

    SCAI_LOG_INFO( logger, "assign done: " << *this );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::buildLocalStorage( _MatrixStorage& storage ) const
{
    if( getColDistribution().isReplicated() )
    {
        // copy local storage with format / value conversion

        storage = *mLocalData;
    }
    else
    {
        // temporary local storage with joined columns needed before

        bool keepDiagonalProperty = true;

        common::shared_ptr<MatrixStorage<ValueType> > tmp( mLocalData->newMatrixStorage() );
        tmp->joinHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), keepDiagonalProperty );
        storage = *tmp;
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::swap( SparseMatrix<ValueType>& other )
{
    Matrix::swapMatrix( other );

    std::swap( mLocalData, other.mLocalData );
    std::swap( mHaloData, other.mHaloData );
    std::swap( mHalo, other.mHalo );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::redistribute( DistributionPtr rowDistributionPtr, DistributionPtr colDistributionPtr )
{
    SCAI_REGION( "Mat.Sp.redistribute" )

    SCAI_LOG_INFO( logger,
                   "redistribute " << *this << ": new row dist = " << *rowDistributionPtr << ", new col dist = " << *colDistributionPtr )

    SCAI_ASSERT_ERROR(
        rowDistributionPtr->getGlobalSize() == mNumRows,
        "size of new row distribution = " << rowDistributionPtr->getGlobalSize() << " does not match number of rows = " << mNumRows );

    SCAI_ASSERT_ERROR(
        colDistributionPtr->getGlobalSize() == mNumColumns,
        *this << ": size of new col distribution = " << colDistributionPtr->getGlobalSize() << " does not match number of columns = " << mNumColumns );

    // Save the current distribution of this matrix; use shared pointers to avoid freeing

    DistributionPtr oldRowDistributionPtr = getRowDistributionPtr();
    DistributionPtr oldColDistributionPtr = getColDistributionPtr();

    // Set the new distributions

    setDistributionPtr( rowDistributionPtr );
    mColDistribution = colDistributionPtr;

    // Handle all cases where we do not have to join the local/halo data of matrix

    if( getRowDistribution() == *oldRowDistributionPtr && getColDistribution() == *oldColDistributionPtr )
    {
        SCAI_LOG_INFO( logger, "row/column distribution are same" )
        return;
    }

    //mLocalData and mHaloData might be recreated so we save theire contextes here to restore them later
    ContextPtr localCtx = mLocalData->getContextPtr();
    ContextPtr haloCtx = mHaloData->getContextPtr();

    if( oldColDistributionPtr->isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, "No column distribution, no halo" )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "remove halo, join local = " << *mLocalData << " and halo = " << *mHaloData )

        bool keepDiagonalProperty = ( *oldColDistributionPtr == *oldRowDistributionPtr );

        mLocalData->joinHalo( *mLocalData, *mHaloData, mHalo, *oldColDistributionPtr, keepDiagonalProperty );

        mHaloData->allocate( mLocalData->getNumRows(), 0 );

        mHalo.clear();

        SCAI_LOG_INFO( logger, "removed column distribution / halo, local data = " << *mLocalData )
    }

    // assign the old local data redistributed to this matrix.

    set( *mLocalData, oldRowDistributionPtr );

    mLocalData->setContextPtr( localCtx );
    mHaloData->setContextPtr( haloCtx );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType> void SparseMatrix<ValueType>::invert( const Matrix& other )
{
    // invert not supported for sparse matrices, so we need a temporary dense matrix

    SCAI_UNSUPPORTED( "invert not supported for sparse matrices, will use a dense matrix" )

    DenseMatrix<ValueType> tmp;

    tmp.invert( other );

    SCAI_LOG_INFO( logger, "tmp = " << tmp << ": inverse of " << other )

    assign( tmp ); // convert back to sparse matrix
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::

set( const MatrixStorage<ValueType>& otherLocalData, DistributionPtr otherDist )
{
    // Note: asssign the other data fitting to the distribution of this matrix

    SCAI_ASSERT_EQUAL_DEBUG( otherDist->getLocalSize(), otherLocalData.getNumRows() )
    SCAI_ASSERT_EQUAL_DEBUG( otherDist->getGlobalSize(), getRowDistribution().getGlobalSize() )

    SCAI_ASSERT_EQUAL_ERROR( otherLocalData.getNumColumns(), getColDistribution().getGlobalSize() )

    if( *otherDist == getRowDistribution() )
    {
        SCAI_LOG_INFO( logger, "same row distribution, assign local" )

        otherLocalData.splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );

    }
    else if( otherDist->isReplicated() )
    {
        // we have to localize and split the other matrix

        SCAI_LOG_INFO( logger, "assign by distribute of replicated data, compute halo" )

        // just split the global replicated data according to
        // column + row distribution of this matrix

        otherLocalData.splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), getRowDistributionPtr().get() );

    }
    else if( getRowDistribution().isReplicated() )
    {
        SCAI_LOG_INFO( logger, "Replication of distributed matrix storage: " << otherLocalData )

        mLocalData->replicate( otherLocalData, *otherDist );

        SCAI_LOG_INFO( logger, "Replicated local storage: " << *mLocalData )

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );

        SCAI_LOG_INFO( logger, "Splitted storage: local = " << *mLocalData << ", halo = " << *mHaloData )
    }
    else
    {
        SCAI_LOG_INFO( logger, "assign is redistribute of distributed matrix" )

        Redistributor redistributor( getRowDistributionPtr(), otherDist );

        SCAI_LOG_INFO( logger,
                       "Redistributor available: source halo = " << redistributor.getHaloSourceSize() << " target halo = " << redistributor.getHaloTargetSize() )

        mLocalData->redistribute( otherLocalData, redistributor );

        SCAI_LOG_INFO( logger, "redistributed, now assign locally" )

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>::~SparseMatrix()
{
    SCAI_LOG_INFO( logger, "~SparseMatrix" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getLocalRow( DenseVector<ValueType>& row, const IndexType iLocal ) const
{
    SCAI_ASSERT_DEBUG( row.getDistribution().isReplicated(), "row vector must be replicated" )

    const Distribution& distributionCol = getColDistribution();

    WriteOnlyAccess<ValueType> rowAccess( row.getLocalValues(), getNumColumns() );

    // Owner of row fills the row by data from local and halo data

    for( IndexType j = 0; j < getNumColumns(); ++j )
    {
        IndexType jLocal = distributionCol.global2local( j );

        SCAI_LOG_TRACE( logger, "global column " << j << " of " << getNumColumns() << " is local " << jLocal )

        if( nIndex != jLocal )
        {
            rowAccess[j] = mLocalData->getValue( iLocal, jLocal );
        }
        else
        {
            // const IndexType jHalo = mHalo.global2halo( j );

            rowAccess[j] = mHaloData->getValue( iLocal, mHalo.global2halo( j ) );
        }

        SCAI_LOG_TRACE( logger, "row[" << j << "] = " << rowAccess[j] )
    }

    // TODO: for optimization make an own loop if distributionCol.isReplicated()
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getDiagonal( Vector& diagonal ) const
{
    if( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( *this << ": set diagonal only supported for row = col distribution." )
    }

    HArray<ValueType> localDiagonal;

    mLocalData->getDiagonal( localDiagonal );

    diagonal.allocate( getRowDistributionPtr() ); // Give the diagonal the right distribution
    diagonal.setValues( localDiagonal ); // Copy values, sizes will fit
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDiagonal( const Vector& diagonal )
{
    if( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( *this << ": set diagonal only supported for row = col distribution." )
    }

    if( getRowDistribution() != diagonal.getDistribution() )
    {
        COMMON_THROWEXCEPTION( diagonal << ": distribution does not fit row distribution of matrix" )
    }

    HArray<ValueType> localDiagonal;

    diagonal.buildValues( localDiagonal );

    // localDiagonal has the same value type as LocalData, no virtual call needed.

    mLocalData->setDiagonalV( localDiagonal );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDiagonal( Scalar value )
{
    if( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    mLocalData->setDiagonal( value.getValue<ValueType>() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::scale( const Vector& scaling )
{
    SCAI_ASSERT_EQUAL( scaling.getDistribution(), getRowDistribution(), "distribution mismatch" )

    HArray<ValueType> localValues;

    scaling.buildValues( localValues );

    mLocalData->scaleRows( localValues );

    // scale Halo storage only if it is used; otherwise there might be a size mismatch

    if( mHaloData->getNumRows() )
    {
        mHaloData->scaleRows( localValues );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::scale( Scalar scaling )
{
    ValueType value = scaling.getValue<ValueType>();

    mLocalData->scale( value );

    if ( mHaloData->getNumRows() )
    {
        mHaloData->scale( value );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::conj()
{
    mLocalData->conj();

    if ( mHaloData->getNumRows() )
    {
        mHaloData->conj();
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesMatrix(
    Matrix& result,
    const Scalar alpha,
    const Matrix& B,
    const Scalar beta,
    const Matrix& C ) const
{
    SCAI_LOG_INFO( logger,
                   "result = alpha * A * B + beta * C with result = " << result << ", alpha = " << alpha << ", A = " << *this << ", B = " << B << ", beta = " << beta << ", C = " << C )

    if( result.getMatrixKind() == Matrix::DENSE )
    {
        // we can deal here with DenseMatrix = SparseMatrix * DenseMatrix + DenseMatrix as it can
        // be considered as matrixTimesVectorN

        DenseMatrix<ValueType>* typedResult = dynamic_cast<DenseMatrix<ValueType>*>( &result );

        SCAI_ASSERT_ERROR( typedResult, "Must be dense matrix<" << getValueType() << "> : " << result )

        const DenseMatrix<ValueType>* typedB = dynamic_cast<const DenseMatrix<ValueType>*>( &B );

        SCAI_ASSERT_ERROR( typedB, "Must be dense matrix<" << getValueType() << "> : " << B )

        ValueType betaVal = beta.getValue<ValueType>();

        const DenseMatrix<ValueType>* typedC = dynamic_cast<const DenseMatrix<ValueType>*>( &C );

        if( betaVal != scai::common::constants::ZERO )
        {
            SCAI_ASSERT_ERROR( typedC, "Must be dense matrix<" << getValueType() << "> : " << C )
        }
        else
        {
            typedC = typedResult; // this alias is well handled
        }

        // Now the typed version can be used

        matrixTimesVectorNImpl( *typedResult, alpha.getValue<ValueType>(), *typedB, beta.getValue<ValueType>(),
                                *typedC );
    }
    else
    {
        SparseMatrix<ValueType>* typedResult = dynamic_cast<SparseMatrix<ValueType>*>( &result );

        SCAI_ASSERT_ERROR( typedResult, "Must be sparse matrix<" << getValueType() << "> : " << result )

        const SparseMatrix<ValueType>* typedB = dynamic_cast<const SparseMatrix<ValueType>*>( &B );

        SCAI_ASSERT_ERROR( typedB, "Must be sparse matrix<" << getValueType() << "> : " << B )

        const SparseMatrix<ValueType>* typedC = dynamic_cast<const SparseMatrix<ValueType>*>( &C );

        SCAI_ASSERT_ERROR( typedC, "Must be sparse matrix<" << getValueType() << "> : " << C )

        // Now the typed version can be used

        typedResult->matrixTimesMatrixImpl( alpha.getValue<ValueType>(), *this, *typedB, beta.getValue<ValueType>(),
                                            *typedC );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixPlusMatrix(
    const Scalar alpha,
    const Matrix& matA,
    const Scalar beta,
    const Matrix& matB )
{
    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", A = " << matA << ", B = " << matB )

    const SparseMatrix<ValueType>* sparseA = dynamic_cast<const SparseMatrix<ValueType>*>( &matA );

    SCAI_ASSERT_ERROR( sparseA, "Must be sparse matrix<" << getValueType() << "> : " << matA )

    const SparseMatrix<ValueType>* sparseB = dynamic_cast<const SparseMatrix<ValueType>*>( &matB );

    SCAI_ASSERT_ERROR( sparseB, "Must be sparse matrix<" << getValueType() << "> : " << matB )

    // Now we can add sparse matrices

    matrixPlusMatrixImpl( alpha.getValue<ValueType>(), *sparseA, beta.getValue<ValueType>(), *sparseB );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixPlusMatrixImpl(
    const ValueType alpha,
    const SparseMatrix<ValueType>& A,
    const ValueType beta,
    const SparseMatrix<ValueType>& B )
{
    SCAI_REGION( "Mat.plusMatrix" )

    // already verified

    SCAI_ASSERT_EQUAL_DEBUG( A.getRowDistribution(), B.getRowDistribution() )
    SCAI_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getColDistribution() )

    if( !B.getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "matrixA * matrixB only supported for replicated columns" << " in matrixB = " << B )
    }

    // Now we can do it completly locally

    Matrix::setDistributedMatrix( A.getRowDistributionPtr(), A.getColDistributionPtr() );

    mLocalData->matrixPlusMatrix( alpha, *A.mLocalData, beta, *B.mLocalData );

    // replicated columns, so no halo needed

    mHaloData->allocate( mNumRows, 0 );
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesMatrixImpl(
    const ValueType alpha,
    const SparseMatrix<ValueType>& A,
    const SparseMatrix<ValueType>& B,
    const ValueType beta,
    const SparseMatrix<ValueType>& C )
{
    SCAI_REGION( "Mat.timesMatrix" )

    SCAI_LOG_DEBUG( logger, "Context lhs before mult " << *(mLocalData->getContextPtr()) )

    if( !B.getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "matrixA * matrixB only supported for replicated columns" << " in matrixB = " << B )
    }

    // already verified

    SCAI_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getRowDistribution() )

    if( beta != scai::common::constants::ZERO )
    {
        SCAI_ASSERT_EQ_ERROR( C.getRowDistribution(), A.getRowDistribution(), "distribution/size mismatch" ) 
        SCAI_ASSERT_EQ_ERROR( C.getColDistribution(), B.getColDistribution(), "distribution/size mismatch" )
    }

    // Now we can do it completly locally

    Matrix::setDistributedMatrix( A.getRowDistributionPtr(), B.getColDistributionPtr() );

    SCAI_LOG_DEBUG( logger, "before matrixTimesMatrix" )
    mLocalData->matrixTimesMatrix( alpha, *A.mLocalData, *B.mLocalData, beta, *C.mLocalData );
    SCAI_LOG_DEBUG( logger, "after matrixTimesMatrix")

    SCAI_LOG_INFO( logger, "local result =  " << *mLocalData )

    if( !A.getColDistribution().isReplicated() )
    {
        CSRStorage<ValueType> haloB; // take CSR format, avoids additional conversions

        // get all needed rows of B, communication plan given by halo schedule of A

        haloB.exchangeHalo( A.getHalo(), B.getLocalStorage(), A.getRowDistribution().getCommunicator() );

        // local = alpha * A_local * B_local + alpha * A_halo * B_halo + C_local

        mLocalData->matrixTimesMatrix( alpha, *A.mHaloData, haloB, static_cast<ValueType>(1.0), *mLocalData );
    }

    // replicated columns, so no halo

    mHaloData->allocate( mNumRows, 0 );
    mHalo.clear();

    SCAI_LOG_DEBUG( logger, "Context lhs after mult " << *(mLocalData->getContextPtr()) )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::haloOperationSync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    HArray<ValueType>& haloX,
    common::function<
    void(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX )> localF,
    common::function<
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX )> haloF ) const
{
    const Communicator& comm = getColDistribution().getCommunicator();

    if( !mHalo.isEmpty() )
    {
        // gather local values of X needed by other processors

        if( mHalo.getProvidesIndexes().size() > 0 )
        {
            // We might receive vaules but do not send them, so the halo might be none empty but provides indexes are.

            SCAI_REGION( "Mat.Sp.syncGatherHalo" )

            SCAI_LOG_INFO( logger,
                           comm << ": gather " << mHalo.getProvidesIndexes().size() << " values of X to provide on " << *localX.getValidContext() );

            HArrayUtils::gather( mTempSendValues, localX, mHalo.getProvidesIndexes() );

            // Note: send values might be fetched to the host by halo exchange
        }

        {
            SCAI_REGION( "Mat.Sp.syncExchangeHalo" )

            SCAI_LOG_INFO( logger, comm << ": halo exchange with : " << mHalo );

            comm.exchangeByPlan( haloX, mHalo.getRequiredPlan(), mTempSendValues, mHalo.getProvidesPlan() );
        }

        if( haloX.size() > 0 )
        {
            SCAI_REGION( "Mat.Sp.syncPrefetchHalo" )

            SCAI_LOG_INFO( logger,
                           comm << ": prefetch " << haloX.size() << " halo values of X to : " << *(mHaloData->getContextPtr()) );

            // During the local computation we prefetch already haloX where it is need

            haloX.prefetch( mHaloData->getContextPtr() );
        }
    }
    else
    {
        SCAI_LOG_INFO( logger, "No halo update needed." )

        haloX.clear(); // sets size to 0, but does not free allocated memory
    }

    {
        SCAI_REGION( "Mat.Sp.syncLocal" )

        SCAI_LOG_INFO( logger,
                       comm << ": synchronous computation localResult[ " << localResult.size() << "] = localF( localMatrix, localX[ " << localX.size() << "] ) on " << *(mLocalData->getContextPtr()) )

        localF( mLocalData.get(), localResult, localX );
    }

    if( haloX.size() > 0 )
    {
        SCAI_REGION( "Mat.Sp.syncHalo" )

        SCAI_LOG_INFO( logger,
                       comm << ": compute with " << haloX.size() << " halo values on " << *(mHaloData->getContextPtr()) );

        // now we can update the result with haloMatrix and halo values of X

        haloF( mHaloData.get(), localResult, haloX );
    }

    SCAI_LOG_DEBUG( logger, "haloOpSync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::vectorHaloOperationSync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    const HArray<ValueType>& localY,
    common::function<
    void(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX )> calcF,
    common::function<
    void(
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX,
        const HArray<ValueType>& localY )> addF ) const
{
    DistributionPtr rowDist = getRowDistributionPtr();
    DistributionPtr colDist = getColDistributionPtr();
    const Communicator& comm = rowDist->getCommunicator();
    IndexType numParts = comm.getSize();
    IndexType myPart = comm.getRank();

    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr localContext = mLocalData->getContextPtr();
    ContextPtr haloContext = mLocalData->getContextPtr();

    IndexType xSize = localX.size();
    SCAI_ASSERT( xSize == rowDist->getLocalSize(),
                 "size mismatch of localX and rowDistribution " << xSize << " != " << rowDist->getLocalSize() )

    IndexType ySize = localY.size();
    IndexType resultSize = localResult.size();
    SCAI_ASSERT( ySize == resultSize, "size mismatch of localY and localResult" << ySize << " != " << resultSize )
    SCAI_ASSERT( ySize == colDist->getLocalSize(),
                 "size mismatch of localY and columnDistribution" << ySize << " != " << colDist->getLocalSize() )

    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;

    // will be done on the host

    std::vector<IndexType> sizes( numParts );
    std::vector<IndexType> offsets;
    comm.allgather( &sizes[0], 1, &xSize );
    offsets = sizes;
    offsets.resize( numParts + 1 );
    sizes2offsets[ hostContext ]( &offsets[0], numParts );

    HArray<ValueType> haloResult( mHalo.getHaloSize() );
    HArray<ValueType> toOthersResult( xSize * numParts );
    HArray<ValueType> fromOthersResult( xSize * numParts );

    if( numParts != 1 )
    {
        // calc halo vector parts
        {
            SCAI_REGION( "Vec.Times.Mat.others" )

            SCAI_LOG_INFO( logger,
                           comm << ": synchronous computation othersResult[ " << toOthersResult.size() << "] = localF( haloMatrix, localX[ " << xSize << "] ) on " << haloContext )

            calcF( mHaloData.get(), haloResult, localX );

            // reassemble halo computation to global indices
            {
                ReadAccess<ValueType> h( haloResult, hostContext );
                WriteAccess<ValueType> o( toOthersResult, hostContext );

                for( IndexType i = 0; i < toOthersResult.size(); ++i )
                {
                    IndexType localIndex = mHalo.global2halo( i );

                    if( localIndex != nIndex )
                    {
                        o[i] = h[localIndex];
                    }
                    else
                    {
                        o[i] = static_cast<ValueType>(0.0);
                    }
                }
            }
        }

        // start vector part exchange
        {
            SCAI_REGION( "vector.swapping" )

            ReadAccess<ValueType> toOthers( toOthersResult, hostContext );
            WriteAccess<ValueType> fromOthers( fromOthersResult, hostContext );

            for( IndexType i = 0; i < numParts; ++i )
            {
                comm.gather( &fromOthers[0], sizes[i], i, &toOthers[offsets[i]] );
            }
        }
        toOthersResult.prefetch( localContext );
    }

    // calc local vector parts
    {
        SCAI_REGION( "Vec.Times.Mat.local" )

        SCAI_LOG_INFO( logger,
                       comm << ": synchronous computation localResult[ " << resultSize << "] = localF( localMatrix, localX[ " << xSize << "] ) on " << localContext )

        calcF( mLocalData.get(), localResult, localX );
    }

    // alpha * ( sum up local vector parts with halo vector parts ) + beta * y
    {
        SCAI_REGION( "Vec.Vec.add" )

        SCAI_LOG_INFO( logger,
                       comm << ": synchronous computation localResult[ " << resultSize << "] = localF( localMatrix, localX[ " << xSize << "] ) on " << localContext )

        if( numParts != 1 )
        {
            ContextPtr contextPtr = Context::getHostPtr();

            WriteAccess<ValueType> localData( localResult, contextPtr );
            ReadAccess<ValueType> otherData( fromOthersResult, contextPtr );

            for( IndexType i = 0; i < numParts; ++i )
            {
                if( i == myPart )
                {
                    continue;
                }

                for( IndexType j = 0; j < xSize; ++j )
                {
                    localData[j] += otherData[i * xSize + j];
                }
            }
        }

        addF( localResult, localResult, localY );
    }

    SCAI_LOG_DEBUG( logger, "vectorHaloOperationSync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::haloOperationAsync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    HArray<ValueType>& haloX,
    common::function<
    tasking::SyncToken*(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX )> localAsyncF,
    common::function<
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX )> haloF ) const
{
    const Communicator& comm = getColDistribution().getCommunicator();

    // We might receive vaules but do not send them, so the halo might be none empty but provides indexes are.

    if( mHalo.getProvidesIndexes().size() > 0 )
    {
        SCAI_REGION( "Mat.Sp.asyncGatherHalo" )

        // gather of halo data cannot be overlapped with local computations on a device
        // Note: gather will be done where denseX is available

        SCAI_LOG_INFO( logger,
                       comm << ": gather " << mHalo.getProvidesIndexes().size() << " values of X to provide on " << *localX.getValidContext() );

        HArrayUtils::gather( mTempSendValues, localX, mHalo.getProvidesIndexes() );

        // prefetch needed otherwise sending will block until local computation has finished

        mTempSendValues.prefetch( comm.getCommunicationContext( mTempSendValues ) );
    }

    common::unique_ptr<tasking::SyncToken> localComputation;

    {
        SCAI_REGION( "Mat.Sp.asyncLocal" )

        SCAI_LOG_INFO( logger,
                       comm << ": start async computation localResult[ " << localResult.size() << "] = localF( localMatrix, localX[ " << localX.size() << "] ) on " << *(mLocalData->getContextPtr()) )

        localComputation.reset( localAsyncF( mLocalData.get(), localResult, localX ) );
    }

    // during local computation we exchange halo and prefetch all needed halo data

    if( !mHalo.isEmpty() )
    {
        SCAI_REGION( "Mat.Sp.asyncExchangeHalo" )

        SCAI_LOG_INFO( logger, comm << ": halo exchange with : " << mHalo );

        comm.exchangeByPlan( haloX, mHalo.getRequiredPlan(), mTempSendValues, mHalo.getProvidesPlan() );

        SCAI_LOG_DEBUG( logger, "Exchange halo done." )
    }

    // start now transfer of the halo values of X to halo context where it is needed

    if( haloX.size() > 0 )

    {
        SCAI_REGION( "Mat.Sp.asyncPrefetchHalo")

        ContextPtr haloLocation = mHaloData->getContextPtr();

        SCAI_LOG_INFO( logger, comm << ": prefetch " << haloX.size() << " halo values of X to : " << *haloLocation );

        haloX.prefetch( haloLocation ); // implicit wait at next access of haloX
    }

    {
        SCAI_REGION( "Mat.Sp.asyncWaitLocal")

        // we must wait for local computation as halo computation updates localResult

        localComputation->wait();

        SCAI_LOG_INFO( logger, comm << ": local computation ready." )
    }

    if( haloX.size() > 0 )
    {
        // now we can update the result with haloMatrix and halo values of X

        SCAI_REGION("Mat.Sp.asyncHalo")

        SCAI_LOG_INFO( logger,
                       comm << ": compute with " << haloX.size() << " halo values on " << *(mHaloData->getContextPtr()) );

        haloF( mHaloData.get(), localResult, haloX );
    }

    SCAI_LOG_DEBUG( logger, "matrixTimesVectorAsync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::vectorHaloOperationAsync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    const HArray<ValueType>& localY,
    common::function<
    tasking::SyncToken*(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX )> calcF,
    common::function<
    /*tasking::SyncToken**/void(
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX,
        const HArray<ValueType>& localY )> addF ) const
{
    DistributionPtr rowDist = getRowDistributionPtr();
    DistributionPtr colDist = getColDistributionPtr();
    const Communicator& comm = rowDist->getCommunicator();
    IndexType numParts = comm.getSize();
    IndexType myPart = comm.getRank();

    ContextPtr hostContext = Context::getHostPtr();
    ContextPtr localContext = mLocalData->getContextPtr();
    ContextPtr haloContext = mLocalData->getContextPtr();

    IndexType xSize = localX.size();
    SCAI_ASSERT( xSize == rowDist->getLocalSize(),
                 "size mismatch of localX and rowDistribution " << xSize << " != " << rowDist->getLocalSize() )

    IndexType ySize = localY.size();
    IndexType resultSize = localResult.size();
    SCAI_ASSERT( ySize == resultSize, "size mismatch of localY and localResult" << ySize << " != " << resultSize )
    SCAI_ASSERT( ySize == colDist->getLocalSize(),
                 "size mismatch of localY and columnDistribution" << ySize << " != " << colDist->getLocalSize() )

    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;

    std::vector<IndexType> sizes( numParts );
    std::vector<IndexType> offsets;
    comm.allgather( &sizes[0], 1, &xSize );
    offsets = sizes;
    offsets.resize( numParts + 1 );
    sizes2offsets[ hostContext ]( &offsets[0], numParts );

    HArray<ValueType> haloResult( mHalo.getHaloSize() );
    HArray<ValueType> toOthersResult( xSize * numParts );
    HArray<ValueType> fromOthersResult( xSize * numParts );

    if( numParts != 1 )
    {
        // calc halo vector parts
        {
            SCAI_REGION( "Vec.Times.Mat.others" )

            SCAI_LOG_INFO( logger,
                           comm << ": asynchronous computation othersResult[" << toOthersResult.size() << "] = localF( haloMatrix, localX[" << xSize << "] ) on " << haloContext )

            delete calcF( mHaloData.get(), haloResult, localX );

            toOthersResult.prefetch( hostContext );

            ReadAccess<ValueType> h( haloResult, hostContext );
            WriteAccess<ValueType> o( toOthersResult, hostContext );

            for( IndexType i = 0; i < toOthersResult.size(); ++i )
            {
                IndexType localIndex = mHalo.global2halo( i );

                if( localIndex != nIndex )
                {
                    o[i] = h[localIndex];
                }
                else
                {
                    o[i] = static_cast<ValueType>(0.0);
                }
            }
        }
    }

    common::unique_ptr<tasking::SyncToken> localComputation;

    // calc local vector parts
    {
        SCAI_REGION( "Vec.Times.Mat.local" )

        SCAI_LOG_INFO( logger,
                       comm << ": asynchronous computation localResult[" << resultSize << "] = localF( localMatrix, localX[" << xSize << "] ) on " << localContext )

        localComputation.reset( calcF( mLocalData.get(), localResult, localX ) );
    }

    ContextPtr contextPtr = Context::getHostPtr();

    if( numParts != 1 )
    {
        // start vector part exchange
        {
            SCAI_REGION( "vector.swapping" )

            SCAI_LOG_INFO( logger,
                           comm << " vector swapping: toOthers[" << toOthersResult.size() << "], fromOthersResult[" << fromOthersResult.size() << "]" )

            ReadAccess<ValueType> toOthers( toOthersResult, contextPtr );
            WriteAccess<ValueType> fromOthers( fromOthersResult, contextPtr );

            for( IndexType i = 0; i < numParts; ++i )
            {
                comm.gather( &fromOthers[0], sizes[i], i, &toOthers[offsets[i]] );
            }
        }

        toOthersResult.prefetch( localContext );
    }

    localComputation->wait();

    // alpha * ( sum up local vector parts with halo vector parts ) + beta * y
    {
        SCAI_REGION( "Vec.Vec.add" )

        SCAI_LOG_INFO( logger,
                       comm << ": asynchronous computation localResult[ " << resultSize << "] = localF( localMatrix, localX[ " << xSize << "] ) on " << localContext )

        if( numParts != 1 )
        {
            WriteAccess<ValueType> localData( localResult, contextPtr );
            ReadAccess<ValueType> otherData( fromOthersResult, contextPtr );

            for( IndexType i = 0; i < numParts; ++i )
            {
                if( i == myPart )
                {
                    continue;
                }

                for( IndexType j = 0; j < xSize; ++j )
                {
                    localData[j] += otherData[i * xSize + j];
                }
            }
        }

        addF( localResult, localResult, localY );
    }

    SCAI_LOG_DEBUG( logger, "vectorHaloOperationAsync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesVectorImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alphaValue,
    const DenseVector<ValueType>& denseX,
    const ValueType betaValue,
    const DenseVector<ValueType>& denseY ) const
{
    SCAI_REGION( "Mat.Sp.timesVector" )

    HArray<ValueType>& localResult = denseResult.getLocalValues();
    const HArray<ValueType>& localY = denseY.getLocalValues();

    const HArray<ValueType>& localX = denseX.getLocalValues();
    HArray<ValueType>& haloX = denseX.getHaloValues();

    // if halo is empty, asynchronous execution is not helpful

    void (scai::lama::MatrixStorage<ValueType>::*matrixTimesVector)(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::matrixTimesVector;

    tasking::SyncToken* (scai::lama::MatrixStorage<ValueType>::*matrixTimesVectorAsync)(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::matrixTimesVectorAsync;

    // routine for halo matrix is same for sync and async version
    using namespace scai::common;

    function<
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX )> haloF =

            // bind( matrixTimesVector, _1, _2, alphaValue, _3, one, cref( localResult ) );

            bind( matrixTimesVector, _1, _2, alphaValue, _3, static_cast<ValueType>(1.0), _2 );

    if( Matrix::SYNCHRONOUS == getCommunicationKind() )
    {
        function<
        void(
            const MatrixStorage<ValueType>* localMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX )> localF =

                bind( matrixTimesVector, _1, _2, alphaValue, _3, betaValue, cref( localY ) );

        haloOperationSync( localResult, localX, haloX, localF, haloF );
    }
    else
    {
        function<
        tasking::SyncToken*(
            const MatrixStorage<ValueType>* localMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX )> localAsyncF =

                bind( matrixTimesVectorAsync, _1, _2, alphaValue, _3, betaValue, cref( localY ) );

        haloOperationAsync( localResult, localX, haloX, localAsyncF, haloF );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::vectorTimesMatrixImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alphaValue,
    const DenseVector<ValueType>& denseX,
    const ValueType betaValue,
    const DenseVector<ValueType>& denseY ) const
{
    SCAI_REGION( "Vec.timesMat.Sp" )

    HArray<ValueType>& localResult = denseResult.getLocalValues();
    const HArray<ValueType>& localX = denseX.getLocalValues();
    const HArray<ValueType>& localY = denseY.getLocalValues();

    // if halo is empty, asynchronous execution is not helpful

    // vectorTimesMatrix: x^ = x * A

    void (MatrixStorage<ValueType>::*vectorTimesMatrix)(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::vectorTimesMatrix;

    tasking::SyncToken* (MatrixStorage<ValueType>::*vectorTimesMatrixAsync)(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::vectorTimesMatrixAsync;

    // vectorPlusVector: result = alpha * x^ + beta * y

    void (*vPlusV)(
        ContextPtr context,
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) = &DenseVector<ValueType>::vectorPlusVector;

    /*tasking::SyncToken**/void (*vPlusVAsync)(
        ContextPtr context,
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) = &DenseVector<ValueType>::vectorPlusVector; //TODO use async: Exception not yet implemented

    // after gather of vector values x^ is on the host
    // todo: think about this if its useful to upload the vector (again)
    ContextPtr hostContext = Context::getHostPtr();

    using namespace scai::common;

    if( Matrix::SYNCHRONOUS == getCommunicationKind() )
    {
        function<
        void(
            const MatrixStorage<ValueType>* localMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX )> calcF = bind( vectorTimesMatrix, _1, _2, static_cast<ValueType>(1.0),
                    _3, static_cast<ValueType>(0.0), _2 );

        function<
        void(
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX,
            const HArray<ValueType>& localY )> addF = bind( vPlusV, hostContext, _1,
                    alphaValue, _2, betaValue, _3 );

        vectorHaloOperationSync( localResult, localX, localY, calcF, addF );
    }
    else
    {
        function<
        tasking::SyncToken*(
            const MatrixStorage<ValueType>* localMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX )> calcAsyncF = bind( vectorTimesMatrixAsync, _1,
                    _2, static_cast<ValueType>(1.0), _3, static_cast<ValueType>(0.0), _2 );

        function<
        /*tasking::SyncToken**/void(
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX,
            const HArray<ValueType>& localY )> addAsyncF = bind( vPlusVAsync, hostContext, _1,
                    alphaValue, _2, betaValue,
                    _3 );

        vectorHaloOperationAsync( localResult, localX, localY, calcAsyncF, addAsyncF );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesVectorNImpl(
    DenseMatrix<ValueType>& result,
    const ValueType alpha,
    const DenseMatrix<ValueType>& x,
    const ValueType beta,
    const DenseMatrix<ValueType>& y ) const
{
    SCAI_REGION( "Mat.Sp.timesVectorN" )

    SCAI_LOG_INFO( logger, "sparseGEMM: " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

    // currently only available for replicated matrices

    SCAI_ASSERT( getRowDistribution().isReplicated(), *this << ": must be replicated" )
    SCAI_ASSERT( getColDistribution().isReplicated(), *this << ": must be replicated" )

    SCAI_ASSERT( x.getRowDistribution().isReplicated(), x << ": must be replicated" )
    SCAI_ASSERT( x.getColDistribution().isReplicated(), x << ": must be replicated" )

    // no alias

    SCAI_ASSERT_DEBUG( &result != &x, "alias of result and X not supported" )

    result.allocate( getRowDistributionPtr(), x.getColDistributionPtr() );

    SCAI_LOG_INFO( logger, "result (allocated) : " << result )

    HArray<ValueType>& resultData = result.getLocalStorage().getData();
    const HArray<ValueType>& xData = x.getLocalStorage().getData();
    const HArray<ValueType>& yData = y.getLocalStorage().getData();

    const IndexType n = result.getLocalStorage().getNumColumns();

    mLocalData->matrixTimesVectorN( resultData, n, alpha, xData, beta, yData );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesVector(
    Vector& result,
    const Scalar alpha,
    const Vector& x,
    const Scalar beta,
    const Vector& y ) const
{
    SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y )

    if( ( &result == &y ) && ( beta != Scalar( 0.0 ) ) )
    {
        SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
    }
    else if( &result == &x )
    {
        COMMON_THROWEXCEPTION( "alias: result = x is not handled, use temporary" )
    }
    else
    {
        // we inherit the row distribution of this matrix to result

        result.allocate( getRowDistributionPtr() );

        // no more to check: result.size() == mNumRows, getRowDistribution() == result.getRowDistribution()
    }

    SCAI_ASSERT_EQUAL_ERROR( x.getDistribution(), getColDistribution() )
    SCAI_ASSERT_EQUAL_ERROR( y.getDistribution(), getRowDistribution() )

    const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
    const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
    DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );

    SCAI_ASSERT( denseX, x << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

    // Note: in case of beta == 0, we might skip this test

    SCAI_ASSERT( denseY, y << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

    SCAI_ASSERT( denseResult, result << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

    matrixTimesVectorImpl( *denseResult, alpha.getValue<ValueType>(), *denseX, beta.getValue<ValueType>(), *denseY );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::vectorTimesMatrix(
    Vector& result,
    const Scalar alpha,
    const Vector& x,
    const Scalar beta,
    const Vector& y ) const
{
    SCAI_LOG_INFO( logger, result << " = " << alpha << " * " << x << " * " << *this << " + " << beta << " * " << y )

    if( ( &result == &y ) && ( beta != Scalar( 0.0 ) ) )
    {
        SCAI_LOG_DEBUG( logger, "alias: result = y is well handled" )
    }
    else if( &result == &x )
    {
        COMMON_THROWEXCEPTION( "alias: result = x is not handled, use temporary" )
    }
    else
    {
        // we inherit the column distribution of this matrix to result

        result.allocate( getColDistributionPtr() );

        // no more to check: result.size() == mNumColumns, getRowDistribution() == result.getColDistribution()
    }

    SCAI_ASSERT_EQUAL_ERROR( x.getDistribution(), getRowDistribution() )
    SCAI_ASSERT_EQUAL_ERROR( y.getDistribution(), getColDistribution() )

    const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
    const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
    DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );

    SCAI_ASSERT( denseX, x << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

    // Note: in case of beta == 0, we might skip this test

    SCAI_ASSERT( denseY, y << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

    SCAI_ASSERT( denseResult, result << ": must be DenseVector<" << common::getScalarType<ValueType>() << ">" )

    vectorTimesMatrixImpl( *denseResult, alpha.getValue<ValueType>(), *denseX, beta.getValue<ValueType>(), *denseY );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesScalar( const Matrix& other, Scalar alpha )
{
    SCAI_LOG_INFO( logger, "this  = " << alpha << " * " << other )

    // should also work fine if other == *this, will not create new data

    assign( other );

    mLocalData->scale( alpha.getValue<ValueType>() );
    if( mHaloData->getNumRows() * mHaloData->getNumColumns() > 0)
    {
		mHaloData->scale( alpha.getValue<ValueType>() );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseMatrix<ValueType>::l1Norm() const
{
    SCAI_REGION( "Mat.Sp.l1Norm" )

    ValueType myValue = mLocalData->l1Norm();
    myValue += mHaloData->l1Norm();


    const Communicator& comm = getRowDistribution().getCommunicator();

    ValueType allValue = comm.sum( myValue );

    SCAI_LOG_INFO( logger, "l1 norm: local value = " << myValue << ", value = " << allValue )

    return Scalar( allValue );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseMatrix<ValueType>::l2Norm() const
{
    SCAI_REGION( "Mat.Sp.l2Norm" )

	ValueType tmp = mLocalData->l2Norm();
    ValueType myValue = tmp * tmp;
    tmp = mHaloData->l2Norm();
	myValue += tmp * tmp;

    const Communicator& comm = getRowDistribution().getCommunicator();

    ValueType allValue = comm.sum( myValue );

	// allValue = ::sqrt( allValue );

	allValue = common::Math::sqrt( allValue );

    SCAI_LOG_INFO( logger, "max norm: local value = " << myValue << ", global value = " << allValue )

    return Scalar( allValue );
}

template<typename ValueType>
Scalar SparseMatrix<ValueType>::maxNorm() const
 {
    SCAI_REGION( "Mat.Sp.maxNorm" )
 
    ValueType myMax = mLocalData->maxNorm();
    ValueType myMaxHalo = mHaloData->maxNorm();
 
    if( myMaxHalo > myMax )
    {
        myMax = myMaxHalo;
    }
 
    const Communicator& comm = getRowDistribution().getCommunicator();
 
    ValueType allMax = comm.max( myMax );

    SCAI_LOG_INFO( logger, "max norm: local max = " << myMax << ", global max = " << allMax )

    return Scalar( allMax );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseMatrix<ValueType>::maxDiffNorm( const Matrix& other ) const
{
    // Implementation works only for same row distribution, replicated col distribution
    // and the same type

    SCAI_REGION( "Mat.Sp.maxDiffNorm" )

    if( ( getRowDistribution() == other.getRowDistribution() ) && getColDistribution().isReplicated()
            && other.getColDistribution().isReplicated() && ( getValueType() == other.getValueType() ) )
    {
        const SparseMatrix<ValueType>* typedOther = dynamic_cast<const SparseMatrix<ValueType>*>( &other );
        SCAI_ASSERT_DEBUG( typedOther, "SERIOUS: wrong dynamic cast: " << other )
        return Scalar( maxDiffNormImpl( *typedOther ) );
    }
    else if( !getColDistribution().isReplicated() )
    {
        // @todo handle maxDiffNorm on sparse matrices with column distribution
        COMMON_THROWEXCEPTION( "maxDiffNorm not available: " << *this << " has column distribution" )
    }
    else
    {
        SCAI_UNSUPPORTED( "maxDiffNorm requires temporary of " << other )
        SparseMatrix<ValueType> typedOther( other, getRowDistributionPtr(), getColDistributionPtr() );
        return Scalar( maxDiffNormImpl( typedOther ) );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseMatrix<ValueType>::maxDiffNormImpl( const SparseMatrix<ValueType>& other ) const
{
    // implementation only supported for same row distributions, replicated columns

    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution(), other.getRowDistribution() )
    SCAI_ASSERT_ERROR( getColDistribution().isReplicated(), *this << ": not replicated column dist" )
    SCAI_ASSERT_ERROR( other.getColDistribution().isReplicated(), other << ": not replicated column dist" )

    ValueType myMaxDiff = mLocalData->maxDiffNorm( other.getLocalStorage() );

    const Communicator& comm = getRowDistribution().getCommunicator();

    ValueType allMaxDiff = comm.max( myMaxDiff );

    SCAI_LOG_INFO( logger, "max diff norm: local max = " << myMaxDiff << ", global max = " << allMaxDiff )

    return allMaxDiff;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseMatrix<ValueType>::getLocalNumValues() const
{
    return mLocalData->getNumValues();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseMatrix<ValueType>::getLocalNumRows() const
{
    return mLocalData->getNumRows();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseMatrix<ValueType>::getLocalNumColumns() const
{
    return mLocalData->getNumColumns();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseMatrix<ValueType>::getNumValues() const
{
    return getRowDistribution().getCommunicator().sum( getPartitialNumValues() );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
IndexType SparseMatrix<ValueType>::getPartitialNumValues() const
{
    return mLocalData->getNumValues() + mHaloData->getNumValues();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseMatrix<ValueType>::getValue( IndexType i, IndexType j ) const
{
    const Distribution& distributionRow = getRowDistribution();
    const Distribution& distributionCol = getColDistribution();
    SCAI_LOG_TRACE( logger, "this(" << i << "," << j << ")" )
    ValueType myValue = static_cast<ValueType>(0.0);
    const IndexType iLocal = distributionRow.global2local( i );

    if( iLocal != nIndex )
    {
        SCAI_LOG_TRACE( logger, "row " << i << " is local " << iLocal )
        IndexType jLocal = distributionCol.global2local( j );

        if( nIndex != jLocal )
        {
            SCAI_LOG_TRACE( logger, "global(" << i << "," << j << ")" " is local(" << iLocal << "," << jLocal << ")" )
            myValue = mLocalData->getValue( iLocal, jLocal );
            SCAI_LOG_TRACE( logger, "found local value " << myValue )
        }
        else
        {
            jLocal = mHalo.global2halo( j );
            SCAI_LOG_TRACE( logger, "global(" << i << "," << j << ")" " is halo(" << iLocal << "," << jLocal << ")" )
            myValue = mHaloData->getValue( iLocal, jLocal );
            SCAI_LOG_TRACE( logger, "found halo value " << myValue )
        }
    }

    SCAI_LOG_TRACE( logger, "myValue = " << myValue )
    myValue = distributionRow.getCommunicator().sum( myValue );
    return Scalar( myValue );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
const Halo& SparseMatrix<ValueType>::getHalo() const
{
    return mHalo;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << getTypeName() << "( size = " << mNumRows << " x " << mNumColumns << ", local = " << *mLocalData
           << ", halo = " << *mHaloData << ", rowdist = " << getRowDistribution() << ", coldist = "
           << getColDistribution() << ")";
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::prefetch() const
{
    mLocalData->prefetch();
    mHaloData->prefetch();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::wait() const
{
    mLocalData->wait();
    mHaloData->wait();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseMatrix<ValueType>::hasDiagonalProperty() const
{
    if( getRowDistribution() != getColDistribution() )
    {
        return false;
    }

    bool localDiagProperty = mLocalData->hasDiagonalProperty();

    bool globalDiagProperty = getRowDistribution().getCommunicator().all( localDiagProperty );

    return globalDiagProperty;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::resetDiagonalProperty()
{
    if( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "diagonal property not possible " )
    }

    this->mLocalData->resetDiagonalProperty();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType SparseMatrix<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
size_t SparseMatrix<ValueType>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>* SparseMatrix<ValueType>::newMatrix() const
{
    COMMON_THROWEXCEPTION( "Can not create a new SparseMatrix with no SparseMatrix format specified" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>* SparseMatrix<ValueType>::copy() const
{
    SCAI_LOG_INFO( logger, "copy of " << *this )

    // copy makes deep copies of the local + halo storage, halo
    //      and uses the same distributions

    shared_ptr<MatrixStorage<ValueType> > newLocalData( mLocalData->copy() );
    shared_ptr<MatrixStorage<ValueType> > newHaloData( mHaloData->copy() );

    SparseMatrix<ValueType>* newSparseMatrix =

        new SparseMatrix<ValueType>( newLocalData, newHaloData, mHalo, getRowDistributionPtr(), getColDistributionPtr() );

    SCAI_LOG_INFO( logger, "copy is " << *newSparseMatrix )

    return newSparseMatrix;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setIdentity( DistributionPtr dist )
{
    allocate( dist, dist );

    const IndexType localNumRows = getRowDistributionPtr()->getLocalSize();

    mLocalData->setIdentity( localNumRows );
    mHaloData->allocate( localNumRows, 0 );

    mHalo.clear(); // no exchange needed

    SCAI_LOG_INFO( logger, *this << ": identity" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDenseData(
    DistributionPtr rowDist,
    DistributionPtr colDist,
    const _HArray& values,
    const Scalar eps )
{
    Matrix::setDistributedMatrix( rowDist, colDist );

    IndexType localNumRows = rowDist->getLocalSize();
    IndexType globalNumCols = colDist->getGlobalSize();

    mLocalData->setDenseData( localNumRows, globalNumCols, values, eps.getValue<ValueType>() );

    if( !colDist->isReplicated() )
    {
        // localize the data according to row distribution, use splitHalo with replicated columns

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
    else
    {
        mHaloData->allocate( localNumRows, 0 );
        mHalo.clear();
    }

    SCAI_LOG_INFO( logger, *this << ": filled by (local) dense data" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setCSRData(
    DistributionPtr rowDist,
    DistributionPtr colDist,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    Matrix::setDistributedMatrix( rowDist, colDist );

    IndexType localNumRows = rowDist->getLocalSize();
    IndexType globalNumCols = colDist->getGlobalSize();

    mLocalData->setCSRData( localNumRows, globalNumCols, numValues, ia, ja, values );

    if( !colDist->isReplicated() )
    {
        // localize the data according to row distribution, use splitHalo with replicated columns

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
    else
    {
        mHaloData->allocate( localNumRows, 0 );
        mHalo.clear();
    }

    SCAI_LOG_INFO( logger, *this << ": filled by (local) dense data" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
size_t SparseMatrix<ValueType>::getMemoryUsage() const
{
    size_t memoryUsage = mLocalData->getMemoryUsage() + mHaloData->getMemoryUsage();

    return getRowDistribution().getCommunicator().sum( memoryUsage );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::writeToFile1(

    const std::string& fileName,
    const File::FileType fileType /* = UNFORMATTED */,
    const common::scalar::ScalarType dataType /* = INTERNAL */,
    const File::IndexDataType indexDataTypeIA /* = LONG */,
    const File::IndexDataType indexDataTypeJA /* = LONG */) const
{
    if( getRowDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file

        const Communicator& comm = getRowDistribution().getCommunicator();

        if( comm.getRank() == 0 )
        {
            mLocalData->writeToFile( fileName, fileType, dataType, indexDataTypeIA, indexDataTypeJA );
        }

        // synchronization to avoid that other processors start with
        // something that might depend on the finally written file

        comm.synchronize();
    }
    else if( hasDiagonalProperty() )
    {
        SCAI_LOG_INFO( logger, "write distributed matrix" )

        const Communicator& comm = getRowDistribution().getCommunicator();

        // as diagonal element is first one we can identify the global id of each row by the column index

        if( getColDistribution().isReplicated() )
        {

            mLocalData->writeToFile( comm.getSize(), comm.getRank(), fileName, fileType, dataType, indexDataTypeIA,
                                     indexDataTypeJA );
        }
        else
        {
            // join the local + halo part before writing it

            CSRStorage<ValueType> local;

            bool keepDiagonalProperty = true;

            local.joinHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), keepDiagonalProperty );

            local.writeToFile( comm.getSize(), comm.getRank(), fileName, fileType, dataType, indexDataTypeIA,
                               indexDataTypeJA );
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( *this << ": write to file not supported with distributions" )
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::readFromFile( const std::string& fileName )
{
    SCAI_REGION( "Mat.Sp.readFromFile" )

    // Take the current default communicator
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    IndexType myRank = comm->getRank();
    IndexType host = 0; // reading processor

    IndexType numRows = 0; // will be the size of the vector
    IndexType numCols = 0; // will be the size of the vector

    if( myRank == host )
    {
        mLocalData->readFromFile( fileName );
        numRows = mLocalData->getNumRows();
        numCols = mLocalData->getNumColumns();
        mHaloData->allocate( numRows, 0 );
        mHalo.clear();
    }

    comm->bcast( &numRows, 1, host );
    comm->bcast( &numCols, 1, host );

    DistributionPtr dist( new CyclicDistribution( numRows, numRows, comm ) );
    DistributionPtr colDist( new NoDistribution( numCols ) );

    if( myRank == host )
    {
        Matrix::setDistributedMatrix( dist, colDist );
    }
    else
    {
        allocate( dist, colDist );
        // other processors have to set sizes of local / halo data
    }
}

template<typename ValueType>
std::string SparseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string("SparseMatrix<") << common::getScalarType<ValueType>() << std::string(">");
    return s.str();
}

template<typename ValueType>
const char* SparseMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return s.c_str();
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( SparseMatrix, ARITHMETIC_HOST_CNT, ARITHMETIC_HOST )

} /* end namespace lama */

} /* end namespace scai */
