/**
 * @file SparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Template specilization of the matrix template for distributed matrixes.
 * @author Jiri Kraus, Thomas Brandes
 * @date 02.04.2012
 */

// hpp
#include <scai/lama/matrix/SparseMatrix.hpp>

// local library

#include <scai/lama/matrix/DenseMatrix.hpp>
#include <scai/lama/matrix/MatrixAssemblyAccess.hpp>
#include <scai/lama/SparseVector.hpp>

#include <scai/lama/storage/MatrixStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>


// scai internal libraries
#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

// std
#include <cmath>
#include <memory>
#include <functional>

using namespace scai::hmemo;
using namespace scai::dmemo;

namespace scai
{

namespace lama
{

using std::shared_ptr;
using std::unique_ptr;
using std::function;
using std::bind;
using std::cref;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using sparsekernel::CSRKernelTrait;

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, SparseMatrix<ValueType>::logger, "Matrix.SparseMatrix" )

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( shared_ptr<MatrixStorage<ValueType> > storage ) :

    Matrix<ValueType>( storage->getNumRows(), storage->getNumColumns() )
{
    mLocalData = storage;

    // create halo ( nLocalRows x 0 ) with same storage format

    mHaloData = shared_ptr<MatrixStorage<ValueType> >( storage->newMatrixStorage() );
    mHaloData->allocate( mLocalData->getNumRows(), 0 );
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( shared_ptr<MatrixStorage<ValueType> > storage, DistributionPtr rowDist ) :

    Matrix<ValueType>(
        rowDist,
        DistributionPtr( new NoDistribution( storage->getNumColumns() ) ) )

{
    mLocalData = storage;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( storage->newMatrixStorage() );

    if ( storage->getNumRows() == rowDist->getLocalSize() )
    {
        // halo is #localRows x 0
        mHaloData->allocate( storage->getNumRows(), 0 );
    }
    else if ( storage->getNumRows() == rowDist->getGlobalSize() )
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
    shared_ptr<MatrixStorage<ValueType> > localData,
    DistributionPtr rowDist,
    DistributionPtr colDist ) :

    Matrix<ValueType>( rowDist, colDist )
{
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getGlobalSize() )
    mLocalData = localData;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( localData->newMatrixStorage() );

    if ( localData->getNumRows() == rowDist->getLocalSize() )
    {
        // build just the halo
        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
    else if ( localData->getNumRows() == rowDist->getGlobalSize() )
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
    shared_ptr<MatrixStorage<ValueType> > localData,
    shared_ptr<MatrixStorage<ValueType> > haloData,
    const Halo& halo,
    DistributionPtr rowDist,
    DistributionPtr colDist ) :

    Matrix<ValueType>( rowDist, colDist )
{
    SCAI_LOG_INFO( logger, "Construct sparse matrix with finalized local, halo storage + Halo" )
    // TODO: asserts for correct sizes of all relevant sizes
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumRows(), rowDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( haloData->getNumRows(), rowDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( haloData->getNumColumns(), halo.getHaloSize() )
    // done by constructor for CRTPMatrix<SparseMatrix<ValueType>, ValueType >:  _Matrix::setDistributedMatrix( rowDist, colDist );
    mLocalData = localData; // flat copy
    mHaloData = haloData; // flat copy
    mHalo = halo;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix() : Matrix<ValueType>( 0, 0 )
{
}

template<typename ValueType>
void SparseMatrix<ValueType>::set(
    shared_ptr<MatrixStorage<ValueType> > localData,
    shared_ptr<MatrixStorage<ValueType> > haloData,
    const Halo& halo,
    DistributionPtr rowDist,
    DistributionPtr colDist ) 
{
    _Matrix::setDistributedMatrix( rowDist, colDist );

    SCAI_LOG_INFO( logger, "Construct sparse matrix with finalized local, halo storage + Halo" )
    // TODO: asserts for correct sizes of all relevant sizes
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumRows(), rowDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( haloData->getNumRows(), rowDist->getLocalSize() )
    SCAI_ASSERT_EQUAL_ERROR( haloData->getNumColumns(), halo.getHaloSize() )

    mLocalData = localData;  // flat copy
    mHaloData = haloData;    // flat copy
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
    _Matrix::checkSettings();

    if ( mHaloData->getNumRows() )
    {
        // If halo storage is available, its size must fit with local storage
        SCAI_ASSERT_EQUAL_DEBUG( mLocalData->getNumRows(), mHaloData->getNumRows() )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
bool SparseMatrix<ValueType>::isConsistent() const
{
    return true;

    const Communicator& comm = getRowDistribution().getCommunicator();

    SCAI_LOG_INFO( logger, comm << ": check for consistency: " << *this )

    IndexType consistencyErrors = 0;

    // ToDo: this implementation should use a corresponding predicate of MatrixStorage

    try
    {
        _Matrix::checkSettings();
        SCAI_ASSERT_EQUAL_ERROR( getRowDistribution().getLocalSize(), mLocalData->getNumRows() )
        SCAI_ASSERT_EQUAL_ERROR( mHaloData->getNumRows(), mLocalData->getNumRows() )
        mLocalData->check( "check for consistency" );
        mHaloData->check( "check for consistency" );
        // ToDo: check Halo
        SCAI_LOG_DEBUG( logger, comm << ": is consistent : " << *this )
    }
    catch ( common::Exception& e )
    {
        SCAI_LOG_INFO( logger, *this << " not consistent: " << e.what() )
        consistencyErrors = 1;
    }

    // use communicator for global reduction to make sure that all processors return same value.

    consistencyErrors = getRowDistribution().getCommunicator().sum( consistencyErrors );

    return 0 == consistencyErrors;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>::SparseMatrix( const SparseMatrix<ValueType>& other ) :

    Matrix<ValueType>( other )

{
    // instead of calling assign( other ), we copy directly to avoid type queries
    // make deep copies of the storage data
    mLocalData = shared_ptr<MatrixStorage<ValueType> >( other.getLocalStorage().copy() );
    mHaloData = shared_ptr<MatrixStorage<ValueType> >( other.getHaloStorage().copy() );
    // just copy the halo
    mHalo = other.getHalo();

    this->setCommunicationKind( other.getCommunicationKind() );
    this->setContextPtr( other.getContextPtr() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::clear()
{
    _Matrix::setReplicatedMatrix( 0, 0 );
    mLocalData->clear();
    mHaloData->clear();
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::purge()
{
    _Matrix::setReplicatedMatrix( 0, 0 );
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
    _Matrix::setDistributedMatrix( distribution, colDistribution );
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
    _Matrix::setReplicatedMatrix( numRows, numColumns );
    mLocalData->allocate( numRows, numColumns );
    mHaloData->allocate( numRows, 0 );
    mHalo.clear();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const _Matrix& matrix )
{
    SCAI_LOG_INFO( logger, "assign " << matrix << " to " << *this )
    this->setContextPtr( matrix.getContextPtr() );
    const SparseMatrix<ValueType>* sparseMatrix = dynamic_cast<const SparseMatrix<ValueType>*>( &matrix );

    if ( sparseMatrix )
    {
        // for a sparse matrix local + halo part can be assigned
        assign( *sparseMatrix );
    }
    else
    {
        // convert dense matrix to a sparse matrix
        DistributionPtr colDist = matrix.getColDistributionPtr();
        DistributionPtr rowDist = matrix.getRowDistributionPtr();
        _Matrix::setDistributedMatrix( rowDist, colDist );
        matrix.buildLocalStorage( *mLocalData ); // local storage with all columns
        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assignTranspose( const _Matrix& matrix )
{
    SCAI_LOG_INFO( logger, "assign transposed " << matrix << " to " << *this )
    const SparseMatrix<ValueType>* sparseMatrix = dynamic_cast<const SparseMatrix<ValueType>*>( &matrix );

    if ( sparseMatrix )
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
    _Matrix::setDistributedMatrix( matrix.getColDistributionPtr(), matrix.getRowDistributionPtr() );

    // Be careful: do not use any more matrix.getRowDistributon or matrix.getColDistribution in case of alias

    if ( getRowDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        SCAI_LOG_DEBUG( logger, "transpose local storage, input = " << matrix.getLocalStorage() )
        mLocalData->assignTranspose( matrix.getLocalStorage() );
        SCAI_LOG_DEBUG( logger, "transposed local storage, is = " << *mLocalData )
        mHaloData->allocate( getRowDistribution().getLocalSize(), 0 );
        mHalo.clear();
        SCAI_LOG_DEBUG( logger, "transposed halo storage, is = " << *mHaloData )
    }
    else if ( getRowDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "transpose not supported for replicated matrices with distributed columns" )
    }
    else if ( getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "transpose not supported for distributed matrices with replicated columns" )
    }
    else
    {
        // rows and columns are both distributed
        SCAI_LOG_INFO( logger, "local transpose of " << matrix.getLocalStorage() )
        mLocalData->assignTranspose( matrix.getLocalStorage() );
        SCAI_LOG_INFO( logger, "local transposed = " << *mLocalData )
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
            const Distribution& dist = getColDistribution();
            const IndexType nJA = sendJA.size();
            WriteAccess<IndexType> ja( sendJA, contextPtr );

            for ( IndexType jj = 0; jj < nJA; ++jj )
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
        mHaloData->setCSRData( mNumLocalRows, getNumColumns(), haloJA.size(), haloRowSizes, haloJA, haloValues );
        // Now build a new halo and localize columns in mHaloData
        mHaloData->buildHalo( mHalo, getColDistribution() );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const SparseMatrix<ValueType>& matrix )
{
    if ( this == &matrix )
    {
        SCAI_LOG_INFO( logger, "self assign sparse matrix = " << matrix )
    }
    else
    {
        SCAI_LOG_INFO( logger, "copy/convert assign sparse matrix = " << matrix )
    }

    _Matrix::setDistributedMatrix( matrix.getRowDistributionPtr(), matrix.getColDistributionPtr() );
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
    _Matrix::setReplicatedMatrix( numRows, numColumns );
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

    if ( !typedStorage )
    {
        mLocalData->assign( storage ); // conversion
        typedStorage = mLocalData.get();
    }

    _Matrix::setDistributedMatrix( rowDist, colDist );
    const Communicator& comm = rowDist->getCommunicator();
    bool isLocal = storage.getNumRows() == rowDist->getLocalSize();
    bool isGlobal = storage.getNumRows() == rowDist->getGlobalSize();
    bool allLocal = comm.all( isLocal ); // local storage on all processors
    bool allGlobal = comm.all( isGlobal ); // global storage on all processors

    if ( allLocal )
    {
        // each processor has exactly its local part
        typedStorage->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
    else if ( allGlobal )
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
    if ( getColDistribution().isReplicated() )
    {
        // copy local storage with format / value conversion
        storage = *mLocalData;
    }
    else
    {
        // temporary local storage with joined columns needed before
        bool keepDiagonalProperty = true;
        shared_ptr<MatrixStorage<ValueType> > tmp( mLocalData->newMatrixStorage() );
        tmp->joinHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), keepDiagonalProperty );
        storage = *tmp;
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::swap( SparseMatrix<ValueType>& other )
{
    _Matrix::swapMatrix( other );
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

    SCAI_ASSERT_EQ_ERROR( rowDistributionPtr->getGlobalSize(), getNumRows(),
                          "size of new row distribution mismatches #rows" );

    SCAI_ASSERT_GE_ERROR( colDistributionPtr->getGlobalSize(), getNumColumns(),
                       "Size of new col distribution = " << *colDistributionPtr << " must be >= #colunns of " << *this );

    // Save the current distribution of this matrix; use shared pointers to avoid freeing

    DistributionPtr oldRowDistributionPtr = getRowDistributionPtr();
    DistributionPtr oldColDistributionPtr = getColDistributionPtr();

    // Set the new distributions

    _Matrix::setDistributedMatrix( rowDistributionPtr, colDistributionPtr );

    // Handle all cases where we do not have to join the local/halo data of matrix

    if ( getRowDistribution() == *oldRowDistributionPtr && getColDistribution() == *oldColDistributionPtr )
    {
        SCAI_LOG_INFO( logger, "row/column distribution are same" )
        return;
    }

    // mLocalData and mHaloData might be recreated so we save their contextes here to restore them later

    ContextPtr localCtx = mLocalData->getContextPtr();
    ContextPtr haloCtx = mHaloData->getContextPtr();

    if ( oldColDistributionPtr->isReplicated() )
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

    // Now we can also increase the number of columns, is tricky here, should be done better

    mLocalData->setDimension( mLocalData->getNumRows(), getNumColumns() );

    // assign the old local data redistributed to this matrix.
    set( *mLocalData, oldRowDistributionPtr );
    mLocalData->setContextPtr( localCtx );
    mHaloData->setContextPtr( haloCtx );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::redistribute( const Redistributor& redistributor, DistributionPtr colDistributionPtr )
{
    SCAI_ASSERT_EQ_ERROR( getRowDistribution(), *redistributor.getSourceDistributionPtr(),
                          "redistributor does not match to actual distribution of this sparse matrix" );

    if ( !getColDistribution().isReplicated() )
    {
        // Halo must be removed before redistribution 

        DistributionPtr repColDistributionPtr( new NoDistribution( getNumColumns() ) );
        redistribute( getRowDistributionPtr(), repColDistributionPtr );
    }

    shared_ptr<MatrixStorage<ValueType> > newData( mLocalData->newMatrixStorage() );
    newData->redistribute( *mLocalData, redistributor );
    mLocalData = newData;

    _Matrix::setDistributedMatrix( redistributor.getTargetDistributionPtr(), getColDistributionPtr() );

    redistribute( getRowDistributionPtr(), colDistributionPtr );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType> void SparseMatrix<ValueType>::invert( const _Matrix& other )
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

    if ( *otherDist == getRowDistribution() )
    {
        SCAI_LOG_INFO( logger, "same row distribution, assign local" )
        otherLocalData.splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
    else if ( otherDist->isReplicated() )
    {
        // we have to localize and split the other matrix
        SCAI_LOG_INFO( logger, "assign by distribute of replicated data, compute halo" )
        // just split the global replicated data according to
        // column + row distribution of this matrix
        otherLocalData.splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), getRowDistributionPtr().get() );
    }
    else if ( getRowDistribution().isReplicated() )
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
                       "Redistributor available: source halo = " << redistributor.getHaloSourceSize() 
                        << " target halo = " << redistributor.getHaloTargetSize() )
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
void SparseMatrix<ValueType>::getLocalRowDense( HArray<ValueType>& row, const IndexType localRowIndex ) const
{
    SCAI_REGION( "Mat.Sp.getLocalRowDense" )

    SCAI_LOG_INFO( logger, "getLocalRowDense( " << localRowIndex << " ) of this matrix: " << *this )

    const Distribution& distributionCol = getColDistribution();

    if ( distributionCol.isReplicated() )
    {
        mLocalData->getRow( row, localRowIndex );
        return;
    }

    row.setSameValue( getNumColumns(), ValueType( 0 ) );

    HArray<ValueType> tmpRow;  // used for row of local, halo data

    // get local part, might be optimized if local part is blocked

    HArray<IndexType> localIndexes;
    distributionCol.getOwnedIndexes( localIndexes );
    mLocalData->getRow( tmpRow, localRowIndex );
    HArrayUtils::scatterImpl( row, localIndexes, true, tmpRow, common::BinaryOp::COPY );

    // get halo part

    mHaloData->getRow( tmpRow, localRowIndex );
    const HArray<IndexType>& haloIndexes = mHalo.getRequiredIndexes();
    HArrayUtils::scatterImpl( row, haloIndexes, true, tmpRow, common::BinaryOp::COPY );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getRowLocal( Vector<ValueType>& row, const IndexType localRowIndex ) const
{
    SCAI_ASSERT_EQ_ERROR( row.getValueType(), getValueType(), "type mismatch" )

    if ( row.getVectorKind() == VectorKind::SPARSE )
    {
        SparseVector<ValueType>& sparseRow = static_cast<SparseVector<ValueType>&>( row );

        sparseRow.allocate( getNumColumns() );

        HArray<ValueType> values;
        HArray<IndexType> indexes;
        getLocalRowSparse( indexes, values, localRowIndex );

        sparseRow.swapSparseValues( indexes, values );
    }
    else if ( row.getVectorKind() == VectorKind::DENSE )
    {
        DenseVector<ValueType>& denseRow = static_cast<DenseVector<ValueType>&>( row );
        denseRow.allocate( getNumColumns() );
        getLocalRowDense( denseRow.getLocalValues(), localRowIndex );
    }
    else 
    {
        COMMON_THROWEXCEPTION( "Unsupported vector kind for row: " << row );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getRow( Vector<ValueType>& row, const IndexType globalRowIndex ) const
{
    SCAI_REGION( "Mat.Sp.getRow" )

    const Distribution& dist = getRowDistribution();
    const Communicator& comm = dist.getCommunicator();

    // if v is not a sparse vector, use a temporary sparse vector

    if ( row.getVectorKind() != VectorKind::SPARSE )
    {
        SCAI_LOG_INFO( logger, "SparseMatrix<" << getValueType() << ">::getRow( DenseVector, " << globalRowIndex << ") requires temporary" )
        SparseVector<ValueType> spRow;
        getRow( spRow, globalRowIndex );
        row.assign( spRow );   // transform the sparse vector into dense vecotr
        return;
    }

    SparseVector<ValueType>& spRow = static_cast<SparseVector<ValueType>&>( row );

    spRow.allocate( getColDistributionPtr() );   // by this way it gets the correct rep distribution

    HArray<IndexType> indexes;   // temporary array for non-zero indexes
    HArray<ValueType> values;    // temporary array for non-zero values

    if ( getRowDistribution().isReplicated() )
    {
        // just pick up my sparse entries from the local part

        mLocalData->getSparseRow( indexes, values, globalRowIndex );
        spRow.swapSparseValues( indexes, values );
        return;
    }

    PartitionId owner = getRowDistribution().findOwner( globalRowIndex );
    
    if ( getColDistribution().isReplicated() )
    {
        SCAI_LOG_INFO( logger, "distributed rows, replicated cols, owner = " << owner << " bcasts row" )

        if ( owner == comm.getRank() )
        {
            // owner gets the sparse row from local data (no halo data) and broadcasts it

            IndexType localRowIndex = getRowDistribution().global2local( globalRowIndex );
            mLocalData->getSparseRow( indexes, values, localRowIndex );
        }

        comm.bcastArray( indexes, owner );  // first call will also bcast the size
        comm.bcastArray( values, indexes.size(), owner );   // 2nd call knows already the size

        spRow.swapSparseValues( indexes, values );

        return;
    }

    // use halo communication 

    if ( owner == comm.getRank() )
    {
        IndexType localRowIndex = getRowDistribution().global2local( globalRowIndex );

        HArray<ValueType> sendData;
        HArray<ValueType> recvData; 

        mHaloData->getRow( sendData, localRowIndex );

        SCAI_LOG_DEBUG( logger, comm << ": as owner send halo row = " << sendData )

        // communicate halo row to other processors corresponding schedule
        // is inverse halo exchange

        const CommunicationPlan recvPlan( NULL, 0 );                  // empty, nothing to receive
        const CommunicationPlan& sendPlan = mHalo.getRequiredPlan();  // only this matters here

        SCAI_LOG_DEBUG( logger, comm << ": owner recvPlan = " << recvPlan << ", sendPlan = " << sendPlan )

        comm.exchangeByPlan( recvData, recvPlan, sendData, sendPlan );

        // copy my local part

        mLocalData->getSparseRow( indexes, values, localRowIndex );

        SCAI_LOG_DEBUG( logger, comm << ": as owner take local row: " << indexes << ", " << values )
    }
    else
    {
        SCAI_LOG_INFO( logger, comm << ": not owner, fill dummy halo, halo size = " << mHaloData->getNumColumns() )

        HArray<ValueType> sendData;   // zero-sized array
        HArray<ValueType> recvData; 

        CommunicationPlan recvPlan;
        recvPlan.extractPlan( mHalo.getProvidesPlan(), owner );
        CommunicationPlan sendPlan( NULL, 0 );                    // empty, nothing to send

        SCAI_LOG_DEBUG( logger, comm << ": not owner recvPlan = " << recvPlan << ", sendPlan = " << sendPlan )

        comm.exchangeByPlan( recvData, recvPlan, sendData, sendPlan );
        
        // only the data received from owner processor is needed

        IndexType n;
        IndexType offset;   // for haloIndexes belonging to owner partition

        mHalo.getProvidesPlan().getInfo( n, offset, owner );

        {
            WriteOnlyAccess<IndexType> wIndexes( indexes, n );
            WriteOnlyAccess<ValueType> wValues( values, n );
            ReadAccess<IndexType> rIndexes( mHalo.getProvidesIndexes() );
            ReadAccess<ValueType> rData( recvData );
          
            IndexType count = 0;  // counts the non-zero values only

            for ( IndexType i = 0; i < n; ++i )
            {
                SCAI_LOG_TRACE( logger, comm << ": got index = " << rIndexes[offset + i] << ", val = " << rData[ i ] )

                if ( rData[i] != common::Constants::ZERO )
                {
                    wIndexes[count] = rIndexes[offset + i];
                    wValues[count] = rData[i];
                    count++;
                }
            }

            wIndexes.resize( count );
            wValues.resize( count );
        }
    }

    spRow.swapSparseValues( indexes, values );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getColumn( Vector<ValueType>& col, const IndexType globalColIndex ) const
{
    SCAI_REGION( "Mat.Sp.getCol" )

    SCAI_LOG_INFO( logger, "SparseMatrix<" << getValueType() << ">::getCol( " << globalColIndex << " )" )

    // if col is not a sparse vector use a temporary sparse vector

    if ( col.getVectorKind() != VectorKind::SPARSE )
    {
        SCAI_LOG_INFO( logger, "SparseMatrix<" << getValueType() << ">::getCol( " << globalColIndex 
                               << ") requires temporary vector for col, kind = " << col.getVectorKind() )
        SparseVector<ValueType> spCol;
        getColumn( spCol, globalColIndex );
        col.assign( spCol );                 // does the required conversion or sparse/dense conversion
        return;
    }

    SCAI_ASSERT_DEBUG( dynamic_cast<SparseVector<ValueType>*>( &col ), "col not SparseVector<" << getValueType() << ">" )

    SparseVector<ValueType>& spCol = static_cast<SparseVector<ValueType>&>( col );

    spCol.allocate( getRowDistributionPtr() );   // resizes local arrays to 0

    // const_cast is safe here as we guarantee consistency with size/distriubiton of spCol

    HArray<IndexType>& rowIndexes = const_cast<HArray<IndexType>&>( spCol.getNonZeroIndexes() );
    HArray<ValueType>& rowValues  = const_cast<HArray<ValueType>&>( spCol.getNonZeroValues() );

    SCAI_ASSERT_EQ_DEBUG( 0, rowIndexes.size(), "allocate of spCol did not clear" );
    SCAI_ASSERT_EQ_DEBUG( 0, rowValues.size(), "allocate of spCol did not clear" );

    IndexType jLocal = getColDistribution().global2local( globalColIndex );

    if ( nIndex != jLocal )
    {
        // column belongs to local storage

        mLocalData->getSparseColumn( rowIndexes, rowValues, jLocal );
    }
    else
    {
        IndexType jHalo = mHalo.global2halo( globalColIndex );

        if ( nIndex != jHalo )
        {
            // column belongs to halo storage

            mHaloData->getSparseColumn( rowIndexes, rowValues, jHalo );
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getLocalRowSparse( HArray<IndexType>& indexes, _HArray& values, const IndexType localRowIndex ) const
{
    SCAI_REGION( "Mat.Sp.getLocalRow" )

    SCAI_LOG_INFO( logger, "getLocalRow( " << localRowIndex << " ) of this matrix: " << *this )

    const Distribution& distributionCol = getColDistribution();

    if ( distributionCol.isReplicated() )
    {
        mLocalData->getSparseRow( indexes, values, localRowIndex );
        return;
    }

    HArray<IndexType> indexes1;
    HArray<IndexType> indexes2;
    HArray<IndexType> haloIndexes2;
    HArray<ValueType> tmpValues;
    HArray<ValueType> values1;
    HArray<ValueType> values2;

    mLocalData->getSparseRow( indexes1, values1, localRowIndex );

    // translate local indexes to global indexes
   
    {
        WriteAccess<IndexType> rIndexes( indexes1 );
        for ( IndexType jj = 0; jj < indexes1.size(); ++jj )
        {
            rIndexes[jj] = distributionCol.local2global( rIndexes[jj] );
        }
    }

    mHaloData->getSparseRow( haloIndexes2, values2, localRowIndex );

    // translate halo indexes to global indexes

    const HArray<IndexType>& haloGlobalIndexes = mHalo.getRequiredIndexes();

    HArrayUtils::gather( indexes2, haloGlobalIndexes, haloIndexes2, common::BinaryOp::COPY );

    ValueType one = 1;  // just union of entries, no scaling
    ValueType zero = 0;  // zero element of sparse arrays

    HArrayUtils::addSparse( indexes, tmpValues, indexes1, values1, zero, one, indexes2, values2, zero, one );

    if ( values.getValueType() == tmpValues.getValueType() )
    {
        HArray<ValueType>& typedValues = static_cast<HArray<ValueType>& >( values );
        typedValues.swap( tmpValues );
    }
    else
    {
        HArrayUtils::assign( values, tmpValues );
    }

    SCAI_LOG_INFO( logger, "getLocalRow( " << localRowIndex << " ) : values = " << values << ", indexes = " << indexes );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setLocalRow( const HArray<ValueType>& row,
        const IndexType localRowIndex,
        const common::BinaryOp op )
{
    SCAI_REGION( "Mat.Sp.setLocalRow" )

    const Distribution& distributionCol = getColDistribution();

    if ( distributionCol.isReplicated() )
    {
        mLocalData->setRow( row, localRowIndex, op );
        return;
    }

    HArray<ValueType> tmpRow;  // used for row of local, halo data

    // set local part

    HArray<IndexType> localIndexes;
    distributionCol.getOwnedIndexes( localIndexes );
    HArrayUtils::gatherImpl( tmpRow, row, localIndexes, common::BinaryOp::COPY );
    mLocalData->setRow( tmpRow, localRowIndex, op );

    // set halo part

    const HArray<IndexType>& haloIndexes = mHalo.getRequiredIndexes();
    HArrayUtils::gatherImpl( tmpRow, row, haloIndexes, common::BinaryOp::COPY );
    mHaloData->setRow( tmpRow, localRowIndex, op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getLocalColumn( HArray<ValueType>& column, const IndexType globalColIndex ) const
{
    SCAI_REGION( "Mat.Sp.getLocalCol" )

    IndexType jLocal = getColDistribution().global2local( globalColIndex );

    const IndexType localRowSize = getRowDistribution().getLocalSize();

    if ( nIndex != jLocal )
    {
        mLocalData->getColumn( column, jLocal );
    }
    else
    {
        IndexType jHalo = mHalo.global2halo( globalColIndex );

        if ( nIndex != jHalo )
        {
            mHaloData->getColumn( column, jHalo );
        }
        else
        {
            column.setSameValue( localRowSize, ValueType( 0 ) );
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setLocalColumn( const HArray<ValueType>& column,
        const IndexType colIndex,
        const common::BinaryOp op )
{
    SCAI_REGION( "Mat.Sp.setLocalCol" )

    const IndexType localRowSize = getRowDistribution().getLocalSize();

    SCAI_ASSERT_EQ_ERROR( column.size(), localRowSize, "serious size mismatch of local column" )

    IndexType jLocal = getColDistribution().global2local( colIndex );

    ReadAccess<ValueType> colAccess( column );

    if ( nIndex != jLocal )
    {
        mLocalData->setColumn( column, jLocal, op );
    }
    else
    {
        IndexType jHalo = mHalo.global2halo( colIndex );

        if ( nIndex != jHalo )
        {
            mHaloData->setColumn( column, jHalo,  op );
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getDiagonal( DenseVector<ValueType>& diagonal ) const
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for square matrices with same row/col distribution" )
    }

    diagonal.allocate( getRowDistributionPtr() );
    mLocalData->getDiagonal( diagonal.getLocalValues() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDiagonal( const DenseVector<ValueType>& diagonal )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "setDiagonal only for square matrices with same row/col distribution" )
    }

    if ( getRowDistribution() != diagonal.getDistribution() )
    {
        COMMON_THROWEXCEPTION( "diagonal must have same distribution as matrix" )
    }

    mLocalData->setDiagonalV( diagonal.getLocalValues() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDiagonal( const ValueType& value )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    mLocalData->setDiagonal( value );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::reduce(
    DenseVector<ValueType>& v, 
    const IndexType dim,
    const common::BinaryOp reduceOp,
    const common::UnaryOp elemOp ) const
{
    SCAI_REGION( "Mat.Sp.reduce" )

    if ( dim == 0 )
    {
        v.allocate( getRowDistributionPtr() );

        // initialize v with neutral element

        v = ValueType( 0 );

        SCAI_LOG_INFO( logger, "reduce dim = 0, v = " << v );

        mLocalData->reduce( v.getLocalValues(), 0, reduceOp, elemOp );
        
        if ( mHaloData->getNumRows() > 0 && mHaloData->getNumColumns() > 0 )
        {
            mHaloData->reduce( v.getLocalValues(), 0, reduceOp, elemOp );
        }

        SCAI_LOG_INFO( logger, "reduced dim = 0, v = " << v );

        return;
    }

    if ( dim == 1 )
    {
        v.allocate( getColDistributionPtr() );

        v = ValueType( 0 );  // neutral element of reduce op

        mLocalData->reduce( v.getLocalValues(), 1, reduceOp, elemOp );

        if ( getRowDistribution().getCommunicator().getSize() == 1 )
        {
            return;   // matrix is replicated, compute just my values
        }

        // rows are distributed

        if ( getColDistribution().getCommunicator().getSize() == 1 )
        {
             SCAI_ASSERT_EQ_ERROR( reduceOp, common::BinaryOp::ADD, "only add supported" )
             getRowDistribution().getCommunicator().sumArray( v.getLocalValues() );
             return;
        }
        
        HArray<ValueType> haloResult;

        if ( !mHalo.isEmpty() )
        {
            const Communicator& comm = getColDistribution().getCommunicator();
            HArray<ValueType> haloData;
            haloData.setSameValue( mHaloData->getNumColumns(), ValueType( 0 ) );
            mHaloData->reduce( haloData, 1, reduceOp, elemOp );
            comm.exchangeByPlan( haloResult, mHalo.getProvidesPlan(), haloData, mHalo.getRequiredPlan() );
        }

        if ( haloResult.size() > 0 )
        {
            HArrayUtils::scatterImpl( v.getLocalValues(), mHalo.getProvidesIndexes(), false, haloResult, reduceOp );
        }

        return;
    }

    COMMON_THROWEXCEPTION( "illegal reduce dim = " << dim )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::scaleRows( const DenseVector<ValueType>& scaleY )
{
    SCAI_ASSERT_EQUAL( scaleY.getDistribution(), getRowDistribution(), "distribution mismatch" )

    const HArray<ValueType>& localValues = scaleY.getLocalValues();

    mLocalData->scaleRows( localValues );

    // scale Halo storage only if it is used; otherwise there might be a size mismatch

    if ( mHaloData->getNumRows() )
    {
        mHaloData->scaleRows( localValues );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::scale( const ValueType& alpha )
{
    mLocalData->scale( alpha );

    if ( mHaloData->getNumRows() )
    {
        mHaloData->scale( alpha );
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
void SparseMatrix<ValueType>::concatenate( 
    dmemo::DistributionPtr rowDist,
    dmemo::DistributionPtr colDist,
    const std::vector<const Matrix<ValueType>*>& matrices )
{
    unique_ptr<_Matrix> mPtr( this->newMatrix() );

    SparseMatrix<ValueType>& newMatrix = static_cast<SparseMatrix<ValueType>& >( *mPtr );
 
    DistributionPtr repColDist = colDist;

    if ( !colDist->isReplicated() )
    {
        repColDist.reset( new NoDistribution( colDist->getGlobalSize() ) );
    }

    newMatrix.allocate( rowDist, repColDist );

    SCAI_LOG_INFO( logger, "Concatenate " << matrices.size() << " matrices into this matrix: " << newMatrix );

    CSRStorage<ValueType> storage;  // reuse in loop

    // Each processor assembles the local part of each input matrix for the result matrix

    {
        MatrixAssemblyAccess<ValueType> assembly( newMatrix );

        IndexType offsetRow = 0;    // row offset where input matrix is set in result matrix
        IndexType offsetCol = 0;    // col offset where input matrix is set in result matrix

        for ( size_t k = 0; k < matrices.size(); ++k )
        {
            const Matrix<ValueType>& m = *matrices[k];

            const Distribution& mRowDist = m.getRowDistribution();

            SCAI_LOG_DEBUG( logger, "dissassemble this matrix: " << m )

            if ( offsetRow + m.getNumRows() > rowDist->getGlobalSize() )
            {
                COMMON_THROWEXCEPTION( "concatenation fails, this arg fails: " << m )
            }

            if ( offsetCol + m.getNumColumns() > colDist->getGlobalSize() )
            {
                COMMON_THROWEXCEPTION( "concatenation fails, arg " << k << " fails: " << m )
            }

            // Build my local part of the current input matrix

            m.buildLocalStorage( storage );

            SCAI_LOG_DEBUG( logger, "Add this storage from matrix " << k << " : " << storage << ", matrix is : " << m )

            ReadAccess<IndexType> rIA( storage.getIA() );
            ReadAccess<IndexType> rJA( storage.getJA() );
            ReadAccess<ValueType> rValues( storage.getValues() );

            // traverse the CSR storage and assemble the entries with the new offsets

            for ( IndexType i = 0; i < storage.getNumRows(); ++i )
            {
                IndexType globalI = mRowDist.local2global( i );

                for ( IndexType jj = rIA[i]; jj < rIA[i+1]; ++jj )
                {
                    IndexType j = rJA[jj];
                    ValueType v = rValues[jj];

                    SCAI_LOG_DEBUG( logger, "push " << offsetRow + globalI << ", " << offsetCol + j << ", " << v )

                    assembly.push( offsetRow + globalI, offsetCol + j, v );
                }
            }

            offsetCol += m.getNumColumns();

            // decide by the sizes where (horizontally or vertically)  we concatenate the next 

            if ( offsetCol == colDist->getGlobalSize() )
            {
                offsetCol = 0;
                offsetRow = offsetRow + m.getNumRows();
            }

            SCAI_LOG_DEBUG( logger, "offsets for next disassembling: row " << offsetRow << ", col " << offsetCol)
        }
    }

    SCAI_LOG_INFO( logger, "assembling finished, new Matrix = " << newMatrix )

    swap( newMatrix );

    redistribute( rowDist, colDist );  // apply column distribution for halo computation
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesMatrix(
    Matrix<ValueType>& result,
    const ValueType alpha,
    const Matrix<ValueType>& B,
    const ValueType beta,
    const Matrix<ValueType>& C ) const
{
    SCAI_LOG_INFO( logger,
                   "result = alpha * A * B + beta * C with result = " << result << ", alpha = " << alpha << ", A = " << *this << ", B = " << B << ", beta = " << beta << ", C = " << C )

    if ( result.getMatrixKind() == MatrixKind::DENSE )
    {
        // we can deal here with DenseMatrix = SparseMatrix * DenseMatrix + DenseMatrix as it can
        // be considered as matrixTimesVectorN

        DenseMatrix<ValueType>& denseResult = static_cast<DenseMatrix<ValueType>&>( result );
 
        SCAI_ASSERT_EQ_ERROR( B.getMatrixKind(), MatrixKind::DENSE, "DenseMatrix = SparseMatrix * B, B must be dense: " << B )

        const DenseMatrix<ValueType>& denseB = static_cast<const DenseMatrix<ValueType>&>( B );

        bool setAlias = false;

        if ( beta != common::Constants::ZERO )
        {
            SCAI_ASSERT_EQ_ERROR( C.getMatrixKind(), MatrixKind::DENSE, "DenseMatrix = SparseMatrix * B + C, C must be dense: " << C )
        }
        else
        {
            // with beta == 0 we do not care at all what type C is

            setAlias = true;
        }

        const DenseMatrix<ValueType>& denseC = setAlias ? denseResult : static_cast<const DenseMatrix<ValueType>&>( C );

        // Now it is just as matrixTimesVector but where Vector is DenseMatrix or just multiple vectors.

        matrixTimesVectorNImpl( denseResult, alpha, denseB, beta, denseC );
    }
    else if ( result.getMatrixKind() == MatrixKind::SPARSE )
    {
        SparseMatrix<ValueType>& sparseResult = static_cast<SparseMatrix<ValueType>&>( result );

        SCAI_ASSERT_EQ_ERROR( B.getMatrixKind(), MatrixKind::SPARSE, "SparseMatrix = SparseMatrix * B, B must be sparse: " << B )
        SCAI_ASSERT_EQ_ERROR( C.getMatrixKind(), MatrixKind::SPARSE, "SparseMatrix = SparseMatrix * SparseMatrix + C, C must be sparse: " << C )

        const SparseMatrix<ValueType>& sparseB = static_cast<const SparseMatrix<ValueType>&>( B );
        const SparseMatrix<ValueType>& sparseC = static_cast<const SparseMatrix<ValueType>&>( C );

        // Now the 'pure' sparse version of matrix mult can be used

        sparseResult.matrixTimesMatrixImpl( alpha, *this, sparseB, beta, sparseC );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Result matrix for matrix-matrix multiplication neither sparse nor dense: " << result )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const Matrix<ValueType>& matA,
    const ValueType beta,
    const Matrix<ValueType>& matB )
{
    SCAI_ASSERT_EQ_ERROR( matA.getRowDistribution(), matB.getRowDistribution(), "size/dist mismatch of matrices to add" )
    SCAI_ASSERT_EQ_ERROR( matB.getColDistribution(), matB.getColDistribution(), "size/dist mismatch of matrices to add" )

    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", A = " << matA << ", B = " << matB )

    SCAI_ASSERT_EQ_ERROR( matA.getMatrixKind(), MatrixKind::SPARSE, "sparseMatrix = alpha * matA + beta * matB, matA must be sparse" )
    SCAI_ASSERT_EQ_ERROR( matB.getMatrixKind(), MatrixKind::SPARSE, "sparseMatrix = alpha * matA + beta * matB, matB must be sparse" )

    const SparseMatrix<ValueType>& sparseA = static_cast<const SparseMatrix<ValueType>&>( matA );
    const SparseMatrix<ValueType>& sparseB = static_cast<const SparseMatrix<ValueType>&>( matB );

    // Now we can add sparse matrices

    matrixPlusMatrixSparse( alpha, sparseA, beta, sparseB );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixPlusMatrixSparse(
    const ValueType alpha,
    const SparseMatrix<ValueType>& A,
    const ValueType beta,
    const SparseMatrix<ValueType>& B )
{
    SCAI_REGION( "Mat.plusMatrix" )

    // already verified

    SCAI_ASSERT_EQ_DEBUG( A.getRowDistribution(), B.getRowDistribution(), "added matrices have different target space" )
    SCAI_ASSERT_EQ_DEBUG( A.getColDistribution(), B.getColDistribution(), "added matrices have different source space" )

    if ( !B.getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "matrixA + matrixB only supported for replicated columns" << " in matrixB = " << B )
    }

    // Now we can do it completly locally
    _Matrix::setDistributedMatrix( A.getRowDistributionPtr(), A.getColDistributionPtr() );
    mLocalData->matrixPlusMatrix( alpha, *A.mLocalData, beta, *B.mLocalData );
    // replicated columns, so no halo needed
    mHaloData->allocate( getNumRows(), 0 );
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
    SCAI_LOG_DEBUG( logger, "Context lhs before mult " << * ( mLocalData->getContextPtr() ) )

    if ( !B.getColDistribution().isReplicated() )
    {
        COMMON_THROWEXCEPTION( "matrixA * matrixB only supported for replicated columns" << " in matrixB = " << B )
    }

    // already verified
    SCAI_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getRowDistribution() )

    if ( beta != common::Constants::ZERO )
    {
        SCAI_ASSERT_EQ_ERROR( C.getRowDistribution(), A.getRowDistribution(), "distribution/size mismatch" )
        SCAI_ASSERT_EQ_ERROR( C.getColDistribution(), B.getColDistribution(), "distribution/size mismatch" )
    }

    // Now we can do it completly locally
    _Matrix::setDistributedMatrix( A.getRowDistributionPtr(), B.getColDistributionPtr() );
    SCAI_LOG_DEBUG( logger, "before matrixTimesMatrix" )
    mLocalData->matrixTimesMatrix( alpha, *A.mLocalData, *B.mLocalData, beta, *C.mLocalData );
    SCAI_LOG_DEBUG( logger, "after matrixTimesMatrix" )
    SCAI_LOG_INFO( logger, "local result =  " << *mLocalData )

    if ( !A.getColDistribution().isReplicated() )
    {
        CSRStorage<ValueType> haloB; // take CSR format, avoids additional conversions
        // get all needed rows of B, communication plan given by halo schedule of A
        haloB.exchangeHalo( A.getHalo(), B.getLocalStorage(), A.getRowDistribution().getCommunicator() );
        // local = alpha * A_local * B_local + alpha * A_halo * B_halo + C_local
        mLocalData->matrixTimesMatrix( alpha, *A.mHaloData, haloB, static_cast<ValueType>( 1.0 ), *mLocalData );
    }

    // replicated columns, so no halo
    mHaloData->allocate( getNumRows(), 0 );
    mHalo.clear();
    SCAI_LOG_DEBUG( logger, "Context lhs after mult " << * ( mLocalData->getContextPtr() ) )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::haloOperationSync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    HArray<ValueType>& haloX,
    function <
    void(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX ) > localF,
    function <
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX ) > haloF ) const
{
    const Communicator& comm = getColDistribution().getCommunicator();

    if ( !mHalo.isEmpty() )
    {
        // gather local values of X needed by other processors
        if ( mHalo.getProvidesIndexes().size() > 0 )
        {
            // We might receive vaules but do not send them, so the halo might be none empty but provides indexes are.
            SCAI_REGION( "Mat.Sp.syncGatherHalo" )
            SCAI_LOG_INFO( logger,
                           comm << ": gather " << mHalo.getProvidesIndexes().size() << " values of X to provide on " << *localX.getValidContext() );
            HArrayUtils::gatherImpl( mTempSendValues, localX, mHalo.getProvidesIndexes(), common::BinaryOp::COPY );
            // Note: send values might be fetched to the host by halo exchange
        }

        {
            SCAI_REGION( "Mat.Sp.syncExchangeHalo" )
            SCAI_LOG_INFO( logger, comm << ": halo exchange with : " << mHalo );
            comm.exchangeByPlan( haloX, mHalo.getRequiredPlan(), mTempSendValues, mHalo.getProvidesPlan() );
        }

        if ( haloX.size() > 0 )
        {
            SCAI_REGION( "Mat.Sp.syncPrefetchHalo" )
            SCAI_LOG_INFO( logger,
                           comm << ": prefetch " << haloX.size() << " halo values of X to : " << * ( mHaloData->getContextPtr() ) );
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
                       comm << ": synchronous computation localResult[ " << localResult.size() << "] = localF( localMatrix, localX[ " << localX.size() << "] ) on " << * ( mLocalData->getContextPtr() ) )
        localF( mLocalData.get(), localResult, localX );
    }

    if ( haloX.size() > 0 )
    {
        SCAI_REGION( "Mat.Sp.syncHalo" )
        SCAI_LOG_INFO( logger,
                       comm << ": compute with " << haloX.size() << " halo values on " << * ( mHaloData->getContextPtr() ) );
        // now we can update the result with haloMatrix and halo values of X
        haloF( mHaloData.get(), localResult, haloX );
    }

    SCAI_LOG_DEBUG( logger, "haloOpSync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::invHaloOperationSync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    HArray<ValueType>& haloX,
    function <
    void(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX ) > localF,
    function <
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX ) > haloF ) const
{
    if ( getRowDistribution().getCommunicator().getSize() == 1 )
    {
        // the full matrix is replicated, the result is distributed, compute just its part

        localF( mLocalData.get(), localResult, localX );

        return;
    }

    const Communicator& comm = getColDistribution().getCommunicator();

    HArray<ValueType> haloResult;  // compute the values for other processors

    if ( !mHalo.isEmpty() )
    {
        SCAI_LOG_INFO( logger, comm << ": compute halo vals for other procs, halo matrix = " << *mHaloData << ", localX = " << localX )

        haloF( mHaloData.get(), haloX, localX );

        SCAI_LOG_DEBUG( logger, comm << ": haloX = " << haloX )

        // send other processors their values, use inverse schedule of halo

        comm.exchangeByPlan( haloResult, mHalo.getProvidesPlan(), haloX, mHalo.getRequiredPlan() );

        SCAI_LOG_DEBUG( logger, comm << ": now exchanged: haloResult = " << haloResult )
    }

    {
        SCAI_REGION( "Mat.Sp.syncLocal" )

        SCAI_LOG_INFO( logger,
                       comm << ": synchronous computation localResult[ " << localResult.size() << "]" <<
                       " = localF( localMatrix, localX[ " << localX.size() << "] )" <<
                       " on " << * ( mLocalData->getContextPtr() ) )

        localF( mLocalData.get(), localResult, localX );
    }

    // Now we have to add the received values from other processors

    if ( haloResult.size() > 0 )
    {
        HArrayUtils::scatterImpl( localResult, mHalo.getProvidesIndexes(), false, haloResult, common::BinaryOp::ADD );
    }

    SCAI_LOG_DEBUG( logger, "invHaloOpSync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::haloOperationAsync(
    HArray<ValueType>& localResult,
    const HArray<ValueType>& localX,
    HArray<ValueType>& haloX,
    function <
    tasking::SyncToken * (
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX ) > localAsyncF,
    function <
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX ) > haloF ) const
{
    const Communicator& comm = getColDistribution().getCommunicator();

    // We might receive vaules but do not send them, so the halo might be none empty but provides indexes are.

    if ( mHalo.getProvidesIndexes().size() > 0 )
    {
        SCAI_REGION( "Mat.Sp.asyncGatherHalo" )
        // gather of halo data cannot be overlapped with local computations on a device
        // Note: gather will be done where denseX is available
        SCAI_LOG_INFO( logger,
                       comm << ": gather " << mHalo.getProvidesIndexes().size() << " values of X to provide on " << *localX.getValidContext() );
        HArrayUtils::gatherImpl( mTempSendValues, localX, mHalo.getProvidesIndexes(), common::BinaryOp::COPY );
        // prefetch needed otherwise sending will block until local computation has finished
        mTempSendValues.prefetch( comm.getCommunicationContext( mTempSendValues ) );
    }

    unique_ptr<tasking::SyncToken> localComputation;
    {
        SCAI_REGION( "Mat.Sp.asyncLocal" )
        SCAI_LOG_INFO( logger,
                       comm << ": start async computation localResult[ " << localResult.size() << "] = localF( localMatrix, localX[ " << localX.size() << "] ) on " << * ( mLocalData->getContextPtr() ) )
        localComputation.reset( localAsyncF( mLocalData.get(), localResult, localX ) );
    }

    // during local computation we exchange halo and prefetch all needed halo data

    if ( !mHalo.isEmpty() )
    {
        SCAI_REGION( "Mat.Sp.asyncExchangeHalo" )
        SCAI_LOG_INFO( logger, comm << ": halo exchange with : " << mHalo );
        comm.exchangeByPlan( haloX, mHalo.getRequiredPlan(), mTempSendValues, mHalo.getProvidesPlan() );
        SCAI_LOG_DEBUG( logger, "Exchange halo done." )
    }

    // start now transfer of the halo values of X to halo context where it is needed

    if ( haloX.size() > 0 )
    {
        SCAI_REGION( "Mat.Sp.asyncPrefetchHalo" )
        ContextPtr haloLocation = mHaloData->getContextPtr();
        SCAI_LOG_INFO( logger, comm << ": prefetch " << haloX.size() << " halo values of X to : " << *haloLocation );
        haloX.prefetch( haloLocation ); // implicit wait at next access of haloX
    }

    {
        SCAI_REGION( "Mat.Sp.asyncWaitLocal" )
        // we must wait for local computation as halo computation updates localResult
        localComputation->wait();
        SCAI_LOG_INFO( logger, comm << ": local computation ready." )
    }

    if ( haloX.size() > 0 )
    {
        // now we can update the result with haloMatrix and halo values of X
        SCAI_REGION( "Mat.Sp.asyncHalo" )
        SCAI_LOG_INFO( logger,
                       comm << ": compute with " << haloX.size() << " halo values on " << * ( mHaloData->getContextPtr() ) );
        haloF( mHaloData.get(), localResult, haloX );
    }

    SCAI_LOG_DEBUG( logger, "matrixTimesVectorAsync done" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesVectorImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alpha,
    const DenseVector<ValueType>& denseX,
    const ValueType beta,
    const DenseVector<ValueType>& denseY ) const
{
    using namespace std;
    using namespace std::placeholders;

    using utilskernel::LArray;

    SCAI_REGION( "Mat.Sp.timesVector" )

    LArray<ValueType>& localResult = denseResult.getLocalValues();
    const LArray<ValueType>& localX = denseX.getLocalValues();
    const LArray<ValueType>& localY = ( beta == common::Constants::ZERO ) ? localResult : denseY.getLocalValues();

    HArray<ValueType>& haloX = denseX.getHaloValues();

    // if halo is empty, asynchronous execution is not helpful

    // matrixTimesVector :  use for bind requires choosing the correct signature here

    void ( scai::lama::MatrixStorage<ValueType>::*matrixTimesVector )(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::matrixTimesVector;

    // matrixTimesVectorAsyn :  use for bind requires choosing the correct signature here

    tasking::SyncToken* ( scai::lama::MatrixStorage<ValueType>::*matrixTimesVectorAsync )(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::matrixTimesVectorAsync;

    // haloF: routine for halo matrix is same for sync and async version

    function <
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX ) > haloF =
            // bind( matrixTimesVector, _1, _2, alpha, _3, one, cref( localResult ) );
            bind( matrixTimesVector, _1, _2, alpha, _3, ValueType( 1 ), _2 );

    if ( SyncKind::SYNCHRONOUS == _Matrix::getCommunicationKind() )
    {
        function <
        void(
            const MatrixStorage<ValueType>* localMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX ) > localF =
                bind( matrixTimesVector, _1, _2, alpha, _3, beta, cref( localY ) );

        haloOperationSync( localResult, localX, haloX, localF, haloF );
    }
    else
    {
        function <
        tasking::SyncToken*(
            const MatrixStorage<ValueType>* localMatrix,
            HArray<ValueType>& localResult,
            const HArray<ValueType>& localX ) > localAsyncF =
                bind( matrixTimesVectorAsync, _1, _2, alpha, _3, beta, cref( localY ) );
        haloOperationAsync( localResult, localX, haloX, localAsyncF, haloF );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::vectorTimesMatrixImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alpha,
    const DenseVector<ValueType>& denseX,
    const ValueType beta,
    const DenseVector<ValueType>& denseY ) const
{
    using namespace std::placeholders;

    SCAI_LOG_INFO( logger, "result = " << alpha << " * x * A + " << beta << " * y"
                   ", x = " << denseX << ", y = " << denseY << ", A = " << *this )

    HArray<ValueType>& localResult = denseResult.getLocalValues();

    const HArray<ValueType>& localX = denseX.getLocalValues();
    const HArray<ValueType>& localY = denseY.getLocalValues();

    void ( MatrixStorage<ValueType>::*vectorTimesMatrix )(
        HArray<ValueType>& result,
        const ValueType alpha,
        const HArray<ValueType>& x,
        const ValueType beta,
        const HArray<ValueType>& y ) const = &MatrixStorage<ValueType>::vectorTimesMatrix;

    using namespace std;

    function <
    void(
        const MatrixStorage<ValueType>* localMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& localX ) > localF =
            bind( vectorTimesMatrix, _1, _2, alpha, _3, beta, cref( localY ) );

    // haloF: localResult = alpha * haloX * haloMatrix

    function <
    void(
        const MatrixStorage<ValueType>* haloMatrix,
        HArray<ValueType>& localResult,
        const HArray<ValueType>& haloX ) > haloF =
            bind( vectorTimesMatrix, _1, _2, alpha, _3, ValueType( 0 ), _2 );

    HArray<ValueType>& haloX = denseX.getHaloValues();  // reuse this array to keep halo values

    invHaloOperationSync( localResult, localX, haloX, localF, haloF );
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

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesScalar( const Matrix<ValueType>& A, ValueType alpha )
{
    SCAI_LOG_INFO( logger, "this  = " << alpha << " * " << A )

    assign( A );  // can also deal with alias

    mLocalData->scale( alpha );

    if ( mHaloData->getNumRows() * mHaloData->getNumColumns() > 0 )
    {
        mHaloData->scale( alpha );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseMatrix<ValueType>::l1Norm() const
{
    SCAI_REGION( "Mat.Sp.l1Norm" )
    RealType<ValueType> myValue = mLocalData->l1Norm();
    myValue += static_cast<RealType<ValueType>>( mHaloData->l1Norm() );
    const Communicator& comm = getRowDistribution().getCommunicator();
    RealType<ValueType> allValue = comm.sum( myValue );
    SCAI_LOG_INFO( logger, "l1 norm: local value = " << myValue << ", value = " << allValue )
    return allValue;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseMatrix<ValueType>::l2Norm() const
{
    SCAI_REGION( "Mat.Sp.l2Norm" )
    RealType<ValueType> tmp = mLocalData->l2Norm();
    RealType<ValueType> myValue = tmp * tmp;
    tmp = mHaloData->l2Norm();
    myValue += tmp * tmp;
    const Communicator& comm = getRowDistribution().getCommunicator();
    RealType<ValueType> allValue = comm.sum( myValue );
    // allValue = ::sqrt( allValue );
    allValue = common::Math::sqrt( allValue );
    SCAI_LOG_INFO( logger, "max norm: local value = " << myValue << ", global value = " << allValue )
    return allValue;
}

template<typename ValueType>
RealType<ValueType> SparseMatrix<ValueType>::maxNorm() const
{
    SCAI_REGION( "Mat.Sp.maxNorm" )

    RealType<ValueType> myMax = mLocalData->maxNorm();
    RealType<ValueType> myMaxHalo = mHaloData->maxNorm();

    if ( myMaxHalo > myMax )
    {
        myMax = myMaxHalo;
    }

    const Communicator& comm = getRowDistribution().getCommunicator();

    RealType<ValueType> allMax = comm.max( myMax );

    SCAI_LOG_INFO( logger, "max norm: local max = " << myMax << ", global max = " << allMax )

    return allMax;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> SparseMatrix<ValueType>::maxDiffNorm( const Matrix<ValueType>& other ) const
{
    // Implementation works only for same row distribution, replicated col distribution
    // and the same type
    SCAI_REGION( "Mat.Sp.maxDiffNorm" )

    bool distributionMatch = getRowDistribution() == other.getRowDistribution();
    bool kindMatch = other.getMatrixKind() == MatrixKind::SPARSE;

    if ( distributionMatch && kindMatch )
    {
        const SparseMatrix<ValueType>* typedOther = dynamic_cast<const SparseMatrix<ValueType>*>( &other );
        SCAI_ASSERT_DEBUG( typedOther, "SERIOUS: wrong dynamic cast: " << other )
        return maxDiffNormImpl( *typedOther );
    }
    else 
    {
        if ( !distributionMatch )
        {
            SCAI_UNSUPPORTED( "maxDiffNorm might be inefficient as distributions do not match" )
        }
        else
        {
            SCAI_UNSUPPORTED( "maxDiffNorm might be inefficient as one matrix is sparse and the other dense" )
        }

        return Matrix<ValueType>::maxDiffNorm( other );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseMatrix<ValueType>::maxDiffNormImpl( const SparseMatrix<ValueType>& other ) const
{
    SCAI_ASSERT_EQ_ERROR( getRowDistribution(), other.getRowDistribution(), "maxDiffNorm: space mismatch" )

    const MatrixStorage<ValueType>* myLocalStorage    = &getLocalStorage();
    const MatrixStorage<ValueType>* otherLocalStorage = &other.getLocalStorage();

    // we might need temporaries for the joined local storage if halo exists

    std::unique_ptr<MatrixStorage<ValueType> > tmp1;
    std::unique_ptr<MatrixStorage<ValueType> > tmp2;

    // ToDo: if col distribution is the same, local and halo might be compared separately.
    //       But halo must be translated to global index space in any case

    if ( !getColDistribution().isReplicated() )
    {
        tmp1.reset( new CSRStorage<ValueType>() );
        tmp1->joinHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), false );
        myLocalStorage = tmp1.get();
    }

    if ( !other.getColDistribution().isReplicated() )
    {
        tmp2.reset( new CSRStorage<ValueType>() );
        tmp2->joinHalo( *other.mLocalData, *other.mHaloData, other.mHalo, other.getColDistribution(), false );
        otherLocalStorage = tmp2.get();
    }

    ValueType myMaxDiff = myLocalStorage->maxDiffNorm( *otherLocalStorage );

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
ValueType SparseMatrix<ValueType>::getValue( IndexType i, IndexType j ) const
{
    const Distribution& distributionRow = getRowDistribution();
    const Distribution& distributionCol = getColDistribution();

    SCAI_LOG_TRACE( logger, "this(" << i << "," << j << ")" )

    ValueType myValue = static_cast<ValueType>( 0 );

    const IndexType iLocal = distributionRow.global2local( i );

    if ( iLocal != nIndex )
    {
        SCAI_LOG_TRACE( logger, "row " << i << " is local " << iLocal )
        IndexType jLocal = distributionCol.global2local( j );

        if ( nIndex != jLocal )
        {
            SCAI_LOG_TRACE( logger, "global(" << i << "," << j << ")" " is local(" << iLocal << "," << jLocal << ")" )
            myValue = mLocalData->getValue( iLocal, jLocal );
            SCAI_LOG_TRACE( logger, "found local value " << myValue )
        }
        else
        {
            jLocal = mHalo.global2halo( j );

            if ( nIndex != jLocal )
            {
                SCAI_LOG_TRACE( logger, "global(" << i << "," << j << ")" " is halo(" << iLocal << "," << jLocal << ")" )
                myValue = mHaloData->getValue( iLocal, jLocal );
                SCAI_LOG_TRACE( logger, "found halo value " << myValue )
            }
        }
    }

    SCAI_LOG_TRACE( logger, "myValue = " << myValue )
    myValue = distributionRow.getCommunicator().sum( myValue );

    return myValue; 
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setValue(
    const IndexType i,
    const IndexType j,
    const ValueType val,
    const common::BinaryOp op )
{
    const Distribution& distributionRow = getRowDistribution();

    const IndexType iLocal = distributionRow.global2local( i );

    if ( iLocal == nIndex )
    {
        return; // this processor does not have the value
    }

    const Distribution& distributionCol = getColDistribution();

    IndexType jLocal = distributionCol.global2local( j );

    if ( nIndex != jLocal )
    {
        mLocalData->setValue( iLocal, jLocal, val, op );
    }
    else
    {
        jLocal = mHalo.global2halo( j );

        if ( nIndex != jLocal )
        {
            mHaloData->setValue( iLocal, jLocal, val, op );
        }
        else
        {
            SCAI_LOG_WARN( logger, "set a non-existing element in sparse matrix ignored" )
        }
    }
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
    stream << getTypeName() << "( size = " << getNumRows() << " x " << getNumColumns() 
           << ", local = " << *mLocalData
           << ", halo = " << *mHaloData;

    // if column and row distribution are equal, write it only once

    if ( getRowDistribution() == getColDistribution() )
    {
        stream << ", dist = " << getRowDistribution() << ")";
    }
    else
    {
        stream << ", rowdist = " << getRowDistribution() 
               << ", coldist = " << getColDistribution() << ")";
    }
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
    if ( getRowDistribution() != getColDistribution() )
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
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "diagonal property not possible " )
    }

    this->mLocalData->resetDiagonalProperty();
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
    _Matrix::setDistributedMatrix( rowDist, colDist );
    IndexType localNumRows = rowDist->getLocalSize();
    IndexType globalNumCols = colDist->getGlobalSize();
    mLocalData->setDenseData( localNumRows, globalNumCols, values, eps.getValue<ValueType>() );

    if ( !colDist->isReplicated() )
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
    _Matrix::setDistributedMatrix( rowDist, colDist );
    IndexType localNumRows = rowDist->getLocalSize();
    IndexType globalNumCols = colDist->getGlobalSize();
    mLocalData->setCSRData( localNumRows, globalNumCols, numValues, ia, ja, values );

    if ( !colDist->isReplicated() )
    {
        // localize the data according to row distribution, use splitHalo with replicated columns
        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
    else
    {
        mHaloData->allocate( localNumRows, 0 );
        mHalo.clear();
    }

    SCAI_LOG_INFO( logger, *this << ": filled by (local) csr data" )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDIAData(
    dmemo::DistributionPtr rowDist,
    dmemo::DistributionPtr colDist,
    const IndexType numDiagonals,
    const hmemo::HArray<IndexType>& offsets,
    const hmemo::_HArray& values )
{
    _Matrix::setDistributedMatrix( rowDist, colDist );
    IndexType localNumRows = rowDist->getLocalSize();
    IndexType globalNumCols = colDist->getGlobalSize();
    mLocalData->setDIAData( localNumRows, globalNumCols, numDiagonals, offsets, values );

    if ( !colDist->isReplicated() )
    {
        // localize the data according to row distribution, use splitHalo with replicated columns
        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
    else
    {
        mHaloData->allocate( localNumRows, 0 );
        mHalo.clear();
    }

    SCAI_LOG_INFO( logger, *this << ": filled by (local) dia data" )
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
std::string SparseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "SparseMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
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

SCAI_COMMON_INST_CLASS( SparseMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
