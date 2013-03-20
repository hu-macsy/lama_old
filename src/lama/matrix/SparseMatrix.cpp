/**
 * @file SparseMatrix.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * $Id$
 */

// hpp
#include <lama/matrix/SparseMatrix.hpp>
#include <lama/matrix/DenseMatrix.hpp>

// others
#include <lama/NoSyncToken.hpp>
#include <lama/LAMAArrayUtils.hpp>

#include <lama/storage/MatrixStorage.hpp>
#include <lama/storage/CSRStorage.hpp>

#include <lama/exception/Exception.hpp>

#include <lama/CommunicatorFactory.hpp>
#include <lama/distribution/NoDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/distribution/Redistributor.hpp>

// tracing
#include <lama/tracing.hpp>

namespace lama
{

using boost::shared_ptr;

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, SparseMatrix<T>::logger, "Matrix.SparseMatrix" );

/* ---------------------------------------------------------------------------------------*/

template<typename T>
SparseMatrix<T>::SparseMatrix( boost::shared_ptr<MatrixStorage<T> > storage ) :

    _SparseMatrix( storage->getNumRows(), storage->getNumColumns() )
{
    mLocalData = storage;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<T> >( storage->create() );
}

/* ---------------------------------------------------------------------------------------*/

template<typename T>
SparseMatrix<T>::SparseMatrix( boost::shared_ptr<MatrixStorage<T> > storage, DistributionPtr rowDist )
    :

    _SparseMatrix( rowDist, DistributionPtr( new NoDistribution( storage->getNumColumns() ) ) )
{
    mLocalData = storage;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<T> >( storage->create() );

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
        LAMA_THROWEXCEPTION( *storage << ": does not fit to row distribution " << *rowDist );
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename T>
SparseMatrix<T>::SparseMatrix(
    boost::shared_ptr<MatrixStorage<T> > localData,
    DistributionPtr rowDist,
    DistributionPtr colDist )
    :

    _SparseMatrix( rowDist, colDist )
{
    LAMA_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getGlobalSize() );

    mLocalData = localData;
    // create empty halo with same storage format
    mHaloData = shared_ptr<MatrixStorage<T> >( localData->create() );

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
        LAMA_THROWEXCEPTION( *localData << ": does not fit to row distribution " << *rowDist );
    }
}

/* ---------------------------------------------------------------------------------------*/

template<typename T>
SparseMatrix<T>::SparseMatrix(
    boost::shared_ptr<MatrixStorage<T> > localData,
    boost::shared_ptr<MatrixStorage<T> > haloData,
    const Halo& halo,
    DistributionPtr rowDist,
    DistributionPtr colDist )
    :

    _SparseMatrix( rowDist, colDist )
{
    LAMA_LOG_INFO( logger, "Construct sparse matrix with finalized local, halo storage + Halo" );

    // TODO: asserts for correct sizes of all relevant sizes

    LAMA_ASSERT_EQUAL_ERROR( localData->getNumRows(), rowDist->getLocalSize() );
    LAMA_ASSERT_EQUAL_ERROR( localData->getNumColumns(), colDist->getLocalSize() );

    LAMA_ASSERT_EQUAL_ERROR( haloData->getNumRows(), rowDist->getLocalSize() );
    LAMA_ASSERT_EQUAL_ERROR( haloData->getNumColumns(), halo.getHaloSize() );

    // done by constructor for _SparseMatrix:  Matrix::setDistributedMatrix( rowDist, colDist );

    mLocalData = localData; // flat copy
    mHaloData = haloData; // flat copy
    mHalo = halo;
}

/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=( const SparseMatrix& matrix )
{
    LAMA_LOG_INFO( logger, " = SparseMatrix : " << matrix );
    assign( matrix );
    return *this;
}

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=( const Matrix& matrix )
{
    LAMA_LOG_INFO( logger, " = Matrix : " << matrix );
    assign( matrix );
    return *this;
}

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=( const Expression<Matrix,Matrix,Times>& exp )
{
    LAMA_LOG_INFO( logger, " = A * B " );
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=( const Expression<Scalar,Matrix,Times>& exp )
{
    LAMA_LOG_INFO( logger, " = alpha * A " );
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=(
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>& exp )
{
    LAMA_LOG_INFO( logger, " = alpha * A * B" );
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=(
    const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> exp )
{
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
SparseMatrix<ValueType>& SparseMatrix<ValueType>::operator=(
    const Expression<Expression<Scalar,Matrix,Times>,Expression<Scalar,Matrix,Times>,Plus> exp )
{
    Matrix::operator=( exp );
    return *this;
}

template<typename ValueType>
void SparseMatrix<ValueType>::checkSettings()
{
    Matrix::checkSettings();

    if( mHaloData->getNumRows() )
    {
        // If halo storage is available, its size must fit with local storage

        LAMA_ASSERT_EQUAL_DEBUG( mLocalData->getNumRows(), mHaloData->getNumRows() );
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

        LAMA_ASSERT_EQUAL_ERROR( getDistribution().getLocalSize(), mLocalData->getNumRows() );
        LAMA_ASSERT_EQUAL_ERROR( mHaloData->getNumRows(), mLocalData->getNumRows() );

        mLocalData->check( "check for consistency" );
        mHaloData->check( "check for consistency" );

        // ToDo: check Halo
    }
    catch( ... )
    {
        consistencyErrors = 1;
    }

    // use communicator for global reduction to make sure that all processors return same value.

    consistencyErrors = getDistribution().getCommunicator().sum( consistencyErrors );

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

    _SparseMatrix( other )

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
void SparseMatrix<ValueType>::allocate( DistributionPtr distribution, DistributionPtr colDistribution )
{
    LAMA_LOG_DEBUG( logger,
                    *this << ": allocate with row dist = " << *distribution << ", col dist = " << *colDistribution );

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
    LAMA_LOG_INFO( logger, "assign " << matrix << " to " << *this );

    this->setContext( matrix.getContextPtr() );

    const _SparseMatrix* sparseMatrix = dynamic_cast<const _SparseMatrix*>( &matrix );

    if( sparseMatrix )
    {
        // for a sparse matrix local + halo part can be assigned

        assign( *sparseMatrix );
    }
    else
    {
        // convert dense matrix to a sparse matrix

        DistributionPtr colDist = matrix.getColDistributionPtr();
        DistributionPtr rowDist = matrix.getDistributionPtr();

        Matrix::setDistributedMatrix( rowDist, colDist );

        matrix.buildLocalStorage( *mLocalData ); // local storage with all columns

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assignTranspose( const Matrix& matrix )
{
    LAMA_LOG_INFO( logger, "assign transposed " << matrix << " to " << *this );

    const SparseMatrix<ValueType>* sparseMatrix = dynamic_cast<const SparseMatrix<ValueType>*>( &matrix );

    if( sparseMatrix )
    {
        assignTransposeImpl( *sparseMatrix );
    }
    else
    {
        LAMA_THROWEXCEPTION( "SparseMatrix::assign currently only implemented for sparse matrices of same type" );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assignTransposeImpl( const SparseMatrix<ValueType>& matrix )

{
    LAMA_LOG_INFO( logger, "transpose sparse matrix with same value type, switch row/col distributions" );

    // assign matrix properties

    Matrix::setDistributedMatrix( matrix.getColDistributionPtr(), matrix.getDistributionPtr() );

    if( getDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        mLocalData->assignTranspose( matrix.getLocalStorage() );
        mHaloData->allocate( getDistribution().getLocalSize(), mNumColumns );
        mHalo.clear();
    }
    else if( getDistribution().isReplicated() )
    {
        LAMA_THROWEXCEPTION( "transpose not supported for replicated matrices with distributed columns" );
    }
    else if( getColDistribution().isReplicated() )
    {
        LAMA_THROWEXCEPTION( "transpose not supported for distributed matrices with replicated columns" );
    }
    else
    {
        // rows and columns are both distributed

        LAMA_LOG_INFO( logger, "local transpose of " << matrix.getLocalStorage() );

        mLocalData->assignTranspose( matrix.getLocalStorage() );

        LAMA_LOG_INFO( logger, "local transposed = " << mLocalData );

        LAMA_LOG_INFO( logger, "halo transpose of " << matrix.getHaloStorage() );

        LAMAArray<IndexType> sendIA;
        LAMAArray<IndexType> sendSizes;
        LAMAArray<IndexType> sendJA;
        LAMAArray<ValueType> sendValues;

        matrix.getHaloStorage().buildCSCData( sendIA, sendJA, sendValues );

        LAMA_LOG_DEBUG( logger,
                        matrix.getHaloStorage() << ": CSC data, IA = " << sendIA << ", JA = " << sendJA << ", Values = " << sendValues );

        // for initial communication we need the sizes and not the offsets

        _MatrixStorage::offsets2sizes( sendSizes, sendIA );

        LAMA_LOG_DEBUG( logger, "sendSizes with " << sendSizes.size() << " entries" );

        LAMAArray<IndexType> recvSizes;

        const Halo& matrixHalo = matrix.getHalo();

        LAMA_ASSERT_EQUAL_DEBUG( sendSizes.size(), matrix.getHaloStorage().getNumColumns() );

        // send all other processors the number of columns

        const Communicator& comm = getDistribution().getCommunicator();

        // Use communication plans of halo in inverse manner to send col sizes

        const CommunicationPlan& sendSizesPlan = matrixHalo.getRequiredPlan();
        const CommunicationPlan& recvSizesPlan = matrixHalo.getProvidesPlan();

        comm.exchangeByPlan( recvSizes, recvSizesPlan, sendSizes, sendSizesPlan );

        // Now we know the sizes, we can pack the data

        LAMAArray<IndexType> recvJA;
        LAMAArray<ValueType> recvValues;

        // Before send of JA: translate back the local (row) indexes to global indexes

        {
            const Distribution& dist = matrix.getDistribution();
            const IndexType nJA = sendJA.size();

            HostWriteAccess<IndexType> ja( sendJA );
            for( IndexType jj = 0; jj < nJA; ++jj )
            {
                ja[jj] = dist.local2global( ja[jj] );
            }
        }

        // Build the communication plan for sparse data

        {
            HostReadAccess<IndexType> sendColSizes( sendSizes );
            HostReadAccess<IndexType> recvColSizes( recvSizes );

            CommunicationPlan sendDataPlan( sendSizesPlan, sendColSizes.get() );
            CommunicationPlan recvDataPlan( recvSizesPlan, recvColSizes.get() );

            comm.exchangeByPlan( recvJA, recvDataPlan, sendJA, sendDataPlan );
            comm.exchangeByPlan( recvValues, recvDataPlan, sendValues, sendDataPlan );
        }

        // Now we have to build the halo by merging the rows of all received data

        LAMAArray<IndexType> haloRowSizes;
        LAMAArray<IndexType> haloJA;
        LAMAArray<ValueType> haloValues;

        const LAMAArray<IndexType>& rowIndexes = matrixHalo.getProvidesIndexes();

        const IndexType mNumLocalRows = getDistribution().getLocalSize();

        {
            // allocate rowSizes as offset array to avoid reallocation
            HostWriteOnlyAccess<IndexType> rowSizes( haloRowSizes, mNumLocalRows + 1 );
        }

        MatrixStorage<ValueType>::joinRows( haloRowSizes, haloJA, haloValues, mNumLocalRows, rowIndexes, recvSizes,
                                            recvJA, recvValues );

        // Now this partition has all sparse data from other partitions to build my new halo

        mHaloData->setCompressThreshold( 1.0 ); // halo rows might be compressed

        _MatrixStorage::sizes2offsets( haloRowSizes );

        mHaloData->setCSRData( mNumLocalRows, mNumColumns, haloJA.size(), haloRowSizes, haloJA, haloValues );

        // Now build a new halo and localize columns in mHaloData

        mHaloData->buildHalo( mHalo, getColDistribution() );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const _SparseMatrix& matrix )
{
    if( this == &matrix )
    {
        LAMA_LOG_INFO( logger, "self assign sparse matrix = " << matrix );
    }
    else
    {
        LAMA_LOG_INFO( logger, "copy/convert assign sparse matrix = " << matrix );
    }

    Matrix::setDistributedMatrix( matrix.getDistributionPtr(), matrix.getColDistributionPtr() );

    // TODO: allow flexibility regarding the context, e.g. format conversion should be done on GPU

    mLocalData->assign( matrix.getLocalStorage() );
    mHaloData->assign( matrix.getHaloStorage() );
    mHalo = matrix.getHalo();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::assign( const _MatrixStorage& storage )
{
    LAMA_LOG_INFO( logger, "assign matrix storage = " << storage );

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
    LAMA_ASSERT_EQUAL_ERROR( storage.getNumColumns(), colDist->getGlobalSize() );

    // split can only be applied to storage of same type, so check if we need to convert value type

    const MatrixStorage<ValueType>* typedStorage = dynamic_cast<const MatrixStorage<ValueType>*>( &storage );

    if( !typedStorage )
    {
        mLocalData->assign( storage ); // conversion
        typedStorage = mLocalData.get();
    }

    Matrix::setDistributedMatrix( rowDist, colDist );

    if( storage.getNumRows() == rowDist->getLocalSize() )
    {
        // build just the halo

        typedStorage->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, NULL );
    }
    else if( storage.getNumRows() == rowDist->getGlobalSize() )
    {
        // we also localize the rows of the matrix

        typedStorage->splitHalo( *mLocalData, *mHaloData, mHalo, *colDist, rowDist.get() );
    }
    else
    {
        LAMA_THROWEXCEPTION( storage << ": does not fit to row distribution " << *rowDist );
    }
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

        boost::shared_ptr<MatrixStorage<ValueType> > tmp( mLocalData->create() );
        tmp->joinHalo( *mLocalData, *mHaloData, mHalo, getColDistribution() );
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
    LAMA_REGION( "Mat.Sp.redistribute" );

    LAMA_LOG_INFO( logger,
                   "redistribute " << *this << ": new row dist = " << *rowDistributionPtr << ", new col dist = " << *colDistributionPtr );

    LAMA_ASSERT_ERROR(
        rowDistributionPtr->getGlobalSize() == mNumRows,
        "size of new row distribution = " << rowDistributionPtr->getGlobalSize() << " does not match number of rows = " << mNumRows );

    LAMA_ASSERT_ERROR(
        colDistributionPtr->getGlobalSize() == mNumColumns,
        *this << ": size of new col distribution = " << colDistributionPtr->getGlobalSize() << " does not match number of columns = " << mNumColumns );

    // Save the current distribution of this matrix; use shared pointers to avoid freeing

    DistributionPtr oldRowDistributionPtr = getDistributionPtr();
    DistributionPtr oldColDistributionPtr = getColDistributionPtr();

    // Set the new distributions

    setDistributionPtr( rowDistributionPtr );
    mColDistribution = colDistributionPtr;

    // Handle all cases where we do not have to join the local/halo data of matrix

    if( getDistribution() == *oldRowDistributionPtr && getColDistribution() == *oldColDistributionPtr )
    {
        LAMA_LOG_INFO( logger, "row/column distribution are same" );
        return;
    }

    //mLocalData and mHaloData might be recreated so we save theire contextes here to restore them later
    ContextPtr localCtx = mLocalData->getContextPtr();
    ContextPtr haloCtx = mHaloData->getContextPtr();

    if( oldColDistributionPtr->isReplicated() )
    {
        LAMA_LOG_DEBUG( logger, "No column distribution, no halo" );
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "remove halo, join local = " << *mLocalData << " and halo = " << *mHaloData );

        mLocalData->joinHalo( *mLocalData, *mHaloData, mHalo, *oldColDistributionPtr );

        mHaloData->allocate( mLocalData->getNumRows(), 0 );

        mHalo.clear();

        LAMA_LOG_INFO( logger, "removed column distribution / halo, local data = " << *mLocalData );
    }

    // assign the old local data redistributed to this matrix.

    set( *mLocalData, oldRowDistributionPtr );

    mLocalData->setContext( localCtx );
    mHaloData->setContext( haloCtx );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType> void SparseMatrix<ValueType>::invert( const Matrix& other )
{
    // invert not supported for sparse matrices, so we need a temporary dense matrix

    LAMA_UNSUPPORTED( "invert not supported for sparse matrices, will use a dense matrix" );

    DenseMatrix<ValueType> tmp;

    tmp.invert( other );

    LAMA_LOG_INFO( logger, "tmp = " << tmp << ": inverse of " << other );

    assign( tmp ); // convert back to sparse matrix
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::

set( const MatrixStorage<ValueType>& otherLocalData, DistributionPtr otherDist )
{
    // Note: asssign the other data fitting to the distribution of this matrix

    LAMA_ASSERT_EQUAL_DEBUG( otherDist->getLocalSize(), otherLocalData.getNumRows() );
    LAMA_ASSERT_EQUAL_DEBUG( otherDist->getGlobalSize(), getDistribution().getGlobalSize() );

    LAMA_ASSERT_EQUAL_ERROR( otherLocalData.getNumColumns(), getColDistribution().getGlobalSize() );

    if( *otherDist == getDistribution() )
    {
        LAMA_LOG_INFO( logger, "same row distribution, assign local" );

        otherLocalData.splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );

    }
    else if( otherDist->isReplicated() )
    {
        // we have to localize and split the other matrix

        LAMA_LOG_INFO( logger, "assign by distribute of replicated data, compute halo" );

        // just split the global replicated data according to
        // column + row distribution of this matrix

        otherLocalData.splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), getDistributionPtr().get() );

    }
    else if( getDistribution().isReplicated() )
    {
        LAMA_LOG_INFO( logger, "Replication of distributed matrix storage: " << otherLocalData );

        mLocalData->replicate( otherLocalData, *otherDist );

        LAMA_LOG_INFO( logger, "Replicated local storage: " << *mLocalData );

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );

        LAMA_LOG_INFO( logger, "Splitted storage: local = " << *mLocalData << ", halo = " << *mHaloData );
    }
    else
    {
        LAMA_LOG_INFO( logger, "assign is redistribute of distributed matrix" );

        Redistributor redistributor( getDistributionPtr(), otherDist );

        LAMA_LOG_INFO( logger,
                       "Redistributor available: source halo = " << redistributor.getHaloSourceSize() << " target halo = " << redistributor.getHaloTargetSize() );

        mLocalData->redistribute( otherLocalData, redistributor );

        LAMA_LOG_INFO( logger, "redistributed, now assign locally" );

        mLocalData->splitHalo( *mLocalData, *mHaloData, mHalo, getColDistribution(), NULL );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SparseMatrix<ValueType>::~SparseMatrix()
{
    LAMA_LOG_INFO( logger, "~SparseMatrix" );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getLocalRow( DenseVector<ValueType>& row, const IndexType iLocal ) const
{
    LAMA_ASSERT_DEBUG( row.getDistribution().isReplicated(), "row vector must be replicated" );

    const Distribution& distributionCol = getColDistribution();

    HostWriteOnlyAccess<ValueType> rowAccess( row.getLocalValues(), getNumColumns() );

    // Owner of row fills the row by data from local and halo data

    for( IndexType j = 0; j < getNumColumns(); ++j )
    {
        IndexType jLocal = distributionCol.global2local( j );

        LAMA_LOG_TRACE( logger, "global column " << j << " of " << getNumColumns() << " is local " << jLocal );

        if( nIndex != jLocal )
        {
            rowAccess[j] = mLocalData->getValue( iLocal, jLocal );
        }
        else
        {
            // const IndexType jHalo = mHalo.global2halo( j );

            rowAccess[j] = mHaloData->getValue( iLocal, mHalo.global2halo( j ) );
        }

        LAMA_LOG_TRACE( logger, "row[" << j << "] = " << rowAccess[j] );
    }

    // TODO: for optimization make an own loop if distributionCol.isReplicated()
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getRow( Vector& row, const IndexType globalRowIndex ) const
{
    LAMA_ASSERT_ERROR( row.getDistribution().isReplicated(), "row vector must be replicated" );

    DenseVector<ValueType>* typedRow = dynamic_cast<DenseVector<ValueType>*>( &row );

    // row must be a DenseVector of same type

    LAMA_ASSERT_ERROR( typedRow, "row is not DenseVector<Matrix::ValueType>" );

    // on a replicated matrix each processor can fill the row

    if( getDistribution().isReplicated() )
    {
        LAMA_LOG_INFO( logger, "get local row " << globalRowIndex );

        getLocalRow( *typedRow, globalRowIndex );
        return;
    }

    // on a distributed matrix, owner fills row and broadcasts it

    const Communicator& comm = getDistribution().getCommunicator();

    // owner fills the row

    IndexType localRowIndex = getDistribution().global2local( globalRowIndex );

    IndexType owner = 0;

    if( localRowIndex != nIndex )
    {
        getLocalRow( *typedRow, localRowIndex );
        owner = comm.getRank() + 1;
        LAMA_LOG_DEBUG( logger,
                        "owner of row " << globalRowIndex << " is " << owner << ", local index = " << localRowIndex );
    }

    owner = comm.sum( owner ) - 1; // get owner via a reduction

    LAMA_ASSERT_ERROR( owner >= 0, "Could not find owner of row " << globalRowIndex );

    {
        HostWriteAccess<ValueType> rowAccess( typedRow->getLocalValues() );
        comm.bcast( rowAccess.get(), getNumColumns(), owner ); // bcast the row
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::getDiagonal( Vector& diagonal ) const
{
    if( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( *this << ": set diagonal only supported for row = col distribution." );
    }

    LAMAArray<ValueType> localDiagonal;

    mLocalData->getDiagonal( localDiagonal );

    diagonal.resize( getDistributionPtr() ); // Give the diagonal the right distribution
    diagonal.setValues( localDiagonal ); // Copy values, sizes will fit
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDiagonal( const Vector& diagonal )
{
    if( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( *this << ": set diagonal only supported for row = col distribution." );
    }

    if( getDistribution() != diagonal.getDistribution() )
    {
        LAMA_THROWEXCEPTION( diagonal << ": distribution does not fit row distribution of matrix" );
    }

    LAMAArray<ValueType> localDiagonal;

    diagonal.buildValues( localDiagonal );

    // localDiagonal has the same value type as LocalData, no virtual call needed.

    mLocalData->setDiagonal( localDiagonal );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setDiagonal( Scalar value )
{
    if( getDistribution() != getColDistribution() ) {
        LAMA_THROWEXCEPTION( "Diagonal calculation only for equal distributions." );
    }

    mLocalData->setDiagonal( value );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::scale( const Vector& scaling )
{
    if( scaling.getDistribution() != getDistribution() )
    {
        LAMA_THROWEXCEPTION(
            scaling << ": scale vector distribution does not fit matrix row distribution " << getDistribution() );
    }

    LAMAArray<ValueType> localValues;

    scaling.buildValues( localValues );

    mLocalData->scale( localValues );

    // scale Halo storage only if it is used; otherwise there might be a size mismatch

    if( mHaloData->getNumRows() )
    {
        mHaloData->scale( localValues );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::scale( Scalar scaling )
{
    if( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "Scale only for equal distributions." );
    }

    mLocalData->scale( scaling );

    if( mHaloData->getNumRows() )
    {
        mHaloData->scale( scaling );
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
    LAMA_LOG_INFO( logger,
                   "result = alpha * A * B + beta * C with result = " << result << ", alpha = " << alpha << ", A = " << *this << ", B = " << B << ", beta = " << beta << ", C = " << C );

    if( result.getMatrixKind() == Matrix::DENSE )
    {
        // we can deal here with DenseMatrix = SparseMatrix * DenseMatrix + DenseMatrix as it can
        // be considered as matrixTimesVectorN

        DenseMatrix<ValueType>* typedResult = dynamic_cast<DenseMatrix<ValueType>*>( &result );

        LAMA_ASSERT_ERROR( typedResult, "Must be dense matrix<" << getValueType() << "> : " << result );

        const DenseMatrix<ValueType>* typedB = dynamic_cast<const DenseMatrix<ValueType>*>( &B );

        LAMA_ASSERT_ERROR( typedB, "Must be dense matrix<" << getValueType() << "> : " << B );

        ValueType betaVal = beta.getValue<ValueType>();

        const DenseMatrix<ValueType>* typedC = dynamic_cast<const DenseMatrix<ValueType>*>( &C );

        if( betaVal != 0 )
        {
            LAMA_ASSERT_ERROR( typedC, "Must be dense matrix<" << getValueType() << "> : " << C );
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

        LAMA_ASSERT_ERROR( typedResult, "Must be sparse matrix<" << getValueType() << "> : " << result );

        const SparseMatrix<ValueType>* typedB = dynamic_cast<const SparseMatrix<ValueType>*>( &B );

        LAMA_ASSERT_ERROR( typedB, "Must be sparse matrix<" << getValueType() << "> : " << B );

        const SparseMatrix<ValueType>* typedC = dynamic_cast<const SparseMatrix<ValueType>*>( &C );

        LAMA_ASSERT_ERROR( typedC, "Must be sparse matrix<" << getValueType() << "> : " << C );

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
    LAMA_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", A = " << matA << ", B = " << matB );

    const SparseMatrix<ValueType>* sparseA = dynamic_cast<const SparseMatrix<ValueType>*>( &matA );

    LAMA_ASSERT_ERROR( sparseA, "Must be sparse matrix<" << getValueType() << "> : " << matA );

    const SparseMatrix<ValueType>* sparseB = dynamic_cast<const SparseMatrix<ValueType>*>( &matB );

    LAMA_ASSERT_ERROR( sparseB, "Must be sparse matrix<" << getValueType() << "> : " << matB );

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
    LAMA_REGION( "Mat.plusMatrix" );

    // already verified

    LAMA_ASSERT_EQUAL_DEBUG( A.getDistribution(), B.getDistribution() );
    LAMA_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getColDistribution() );

    if( !B.getColDistribution().isReplicated() )
    {
        LAMA_THROWEXCEPTION( "matrixA * matrixB only supported for replicated columns" << " in matrixB = " << B );
    }

    // Now we can do it completly locally

    Matrix::setDistributedMatrix( A.getDistributionPtr(), A.getColDistributionPtr() );

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
    LAMA_REGION( "Mat.timesMatrix" );

    LAMA_LOG_DEBUG( logger, "Context lhs before mult " << mLocalData->getContext() );

    if( !B.getColDistribution().isReplicated() )
    {
        LAMA_THROWEXCEPTION( "matrixA * matrixB only supported for replicated columns" << " in matrixB = " << B );
    }

    // already verified

    LAMA_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getDistribution() );

    if( beta != 0.0 )
    {
        LAMA_ASSERT_ERROR( C.getDistribution() == A.getDistribution(),
                           "Row distribution must be " << A.getDistribution() << ": " << C );
        LAMA_ASSERT_ERROR( C.getColDistribution() == A.getColDistribution(),
                           "Col distribution must be " << A.getColDistribution() << ": " << C );
    }

    // Now we can do it completly locally

    Matrix::setDistributedMatrix( A.getDistributionPtr(), B.getColDistributionPtr() );

    mLocalData->matrixTimesMatrix( alpha, *A.mLocalData, *B.mLocalData, beta, *C.mLocalData );

    LAMA_LOG_INFO( logger, "local result =  " << *mLocalData );

    if( !A.getColDistribution().isReplicated() )
    {
        CSRStorage<ValueType> haloB; // take CSR format, avoids additional conversions

        // get all needed rows of B, communication plan given by halo schedule of A

        haloB.exchangeHalo( A.getHalo(), B.getLocalStorage(), A.getDistribution().getCommunicator() );

        // local = alpha * A_local * B_local + alpha * A_halo * B_halo + C_local

        mLocalData->matrixTimesMatrix( alpha, *A.mHaloData, haloB, 1.0, *mLocalData );
    }

    // replicated columns, so no halo

    mHaloData->allocate( mNumRows, 0 );
    mHalo.clear();

    LAMA_LOG_DEBUG( logger, "Context lhs after mult " << mLocalData->getContext() );
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
    LAMA_REGION( "Mat.Sp.timesVectorImpl" );

    ContextPtr localContext = mLocalData->getContextPtr();

    //Prefetch matrix to local location while exchanging the halo
    //see comment below why prefetch for halo location is not started here
    LAMA_LOG_DEBUG( logger, "Starting prefetch of input data for local computation to: " << localContext );

    mLocalData->prefetch();

    //It makes no sense to prefetch denseX because, if a transfer is started
    //the halo update needs to wait for this transfer to finish
    if( betaValue != zero )
    {
        denseY.prefetch( localContext );
    }

    if( SYNCHRONOUS == getCommunicationKind() && !mHalo.isEmpty() )
    {
        LAMA_REGION( "Mat.Sp.timesVector::syncUpdateHalo");
        //1. gather
        // We might receive vaules but do not send them, so the halo might be none empty but provides indexes are.
        if( mHalo.getProvidesIndexes().size() > 0 )
        {
            LAMAArrayUtils::gather( mTempSendValues, denseX.getLocalValues(), mHalo.getProvidesIndexes() );
        }
        //2. exchange by plan
        getColDistribution().getCommunicator().exchangeByPlan( denseX.getHaloValues(), mHalo.getRequiredPlan(),
                mTempSendValues, mHalo.getProvidesPlan() );
        //denseX.updateHalo( mHalo );
        LAMA_LOG_INFO( logger, "Synchronous update of halo values done." );
    }
    else
    {
        LAMA_LOG_INFO( logger, "No update for halo values necessary." );
    }

    const LAMAArray<ValueType>& localX = denseX.getLocalValues();
    const LAMAArray<ValueType>& localY = denseY.getLocalValues();
    LAMAArray<ValueType>& localResult = denseResult.getLocalValues();

    LAMA_LOG_INFO( logger, "Starting local computation." );

    //get best possible routine for computation of local matrix

    LAMA_LOG_DEBUG( logger,
                    " Calling LAMAInterface for z = alpha * A * x + beta * y" << " with A = " << mLocalData << ", x = " << localX << ", y = " << localY << ", z = " << localResult );

    std::auto_ptr<SyncToken> localComputation;

    if( ASYNCHRONOUS == getCommunicationKind() )
    {
        // We might receive vaules but do not send them, so the halo might be none empty but provides indexes are.
        if( mHalo.getProvidesIndexes().size() > 0 )
        {
            //1. gather
            LAMAArrayUtils::gather( mTempSendValues, denseX.getLocalValues(), mHalo.getProvidesIndexes() );
            //localX.prefetch( ContextFactory::getContext( Context::Host ) );
            //prefetch neede because if the copy is started after the computation the copy blocks
            mTempSendValues.prefetch( ContextFactory::getContext( Context::Host ) );
        }
        LAMA_REGION( "Mat.Sp.timesVector::localAsync" );
        LAMA_LOG_INFO( logger, "Starting asynchronous computation of local values on " << *localContext );
        localComputation = mLocalData->matrixTimesVectorAsync( localResult, alphaValue, localX, betaValue, localY );
    }
    else
    {
        LAMA_REGION( "Mat.Sp.timesVector::localSync" );
        LAMA_LOG_INFO( logger, "Starting synchronous computation of local values on " << *localContext );
        mLocalData->matrixTimesVector( localResult, alphaValue, localX, betaValue, localY );
        localComputation.reset( new NoSyncToken() );
    }

    // prefetch matrix to halo location while waiting for the halo exchange
    // this prefetch is not started together with the local location prefetch because only on running
    // prefetch is possible at any time. So
    // pefetch( localContext, false );
    // prefetch( mHaloLocation, false );
    // would lead to an immeadiate synchronize for the local location prefetch

    mHaloData->prefetch();

    if( ASYNCHRONOUS == getCommunicationKind() && !mHalo.isEmpty() )
    {
        LAMA_REGION( "Mat.Sp.timesVector::updateHalo");
        //2. do exchange by plan
        getColDistribution().getCommunicator().exchangeByPlan( denseX.getHaloValues(), mHalo.getRequiredPlan(),
                mTempSendValues, mHalo.getProvidesPlan() );
        //denseX.updateHalo( mHalo );
        LAMA_LOG_INFO( logger, "Synchronous update of halo values parallel to asynchronous local computation done." );
    }

    LAMA_LOG_INFO( logger, "Halo available." );

    if( mHalo.getHaloSize() > 0 )
    {
        /*  @todo: this check is only useful for CSR and ELL, not for DIA, COO, JDS
         @todo: even for CSR and ELL the Halo is probably not always sparse

         if ( mHaloData->getRowIndexes().size() == 0 )
         {
         LAMA_LOG_FATAL( logger, "PERFORMANCE WARNING: mHaloData->mRowIndexes is not correctly initialized."
         << " mHaloData = " << *mHaloData << " mHaloData->mRowIndexes.size() ==  "
         << mHaloData->getRowIndexes().size() );
         }
         */

        LAMA_REGION( "Mat.Sp.timesVector::computeHalo" );

        LAMA_LOG_INFO( logger, "Halo compute, size = " << mHalo.getHaloSize() );

        const LAMAArray<ValueType>& haloX = denseX.getHaloValues();

        ContextPtr haloLocation = mHaloData->getContextPtr();

        haloX.prefetch( haloLocation );

        LAMA_LOG_INFO( logger, "Halo compute, haloX = " << haloX );

        LAMA_LOG_DEBUG( logger,
                        " Calling LAMAInterface for z += alpha * A * x + beta * y " << " with A = " << mHaloData << ", x = " << haloX << ", y = " << localY << ", z = " << localResult );

        LAMA_LOG_INFO( logger, "Starting halo computation after the local computation on " << *haloLocation );
        {
            LAMA_REGION( "Mat.Sp.timesVector::waitLocal");
            localComputation->wait();
        }
        LAMA_LOG_INFO( logger, "Local computation done." );
        {
            LAMA_REGION("Mat.Sp.timesVector::halo_timesVector");
            mHaloData->matrixTimesVector( localResult, alphaValue, haloX, 1.0, localResult );
        }
    }
    {
        LAMA_REGION( "Mat.Sp.timesVector::waitLocalNoHalo");
        // Need to synchronize the localComputation if halo is empty and we have async comp
        localComputation->wait();
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
    LAMA_REGION( "Mat.Sp.timesVectorN" );

    LAMA_LOG_INFO( logger, "sparseGEMM: " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y );

    // currently only available for replicated matrices

    LAMA_ASSERT( getDistribution().isReplicated(), *this << ": must be replicated" );
    LAMA_ASSERT( getColDistribution().isReplicated(), *this << ": must be replicated" );

    LAMA_ASSERT( x.getDistribution().isReplicated(), x << ": must be replicated" );
    LAMA_ASSERT( x.getColDistribution().isReplicated(), x << ": must be replicated" );

    // no alias

    LAMA_ASSERT_DEBUG( &result != &x, "alias of result and X not supported" );

    result.allocate( getDistributionPtr(), x.getColDistributionPtr() );

    LAMA_LOG_INFO( logger, "result (allocated) : " << result );

    LAMAArray<ValueType>& resultData = result.getLocalStorage().getData();
    const LAMAArray<ValueType>& xData = x.getLocalStorage().getData();
    const LAMAArray<ValueType>& yData = y.getLocalStorage().getData();

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
    LAMA_REGION( "Mat.Sp.timesVector" );

    LAMA_LOG_INFO( logger, result << " = " << alpha << " * " << *this << " * " << x << " + " << beta << " * " << y );

    if( ( &result == &y ) && ( beta != Scalar( 0.0 ) ) )
    {
        LAMA_LOG_DEBUG( logger, "alias: result = y is well handled" );
    }
    else if( &result == &x )
    {
        LAMA_THROWEXCEPTION( "alias: result = x is not handled, use temporary" );
    }
    else
    {
        // we inherit the row distribution of this matrix to result

        result.resize( getDistributionPtr() );

        // no more to check: result.size() == mNumRows, getDistribution() == result.getDistribution()
    }

    LAMA_ASSERT_EQUAL( x.getDistribution(), getColDistribution() );
    LAMA_ASSERT_EQUAL( y.getDistribution(), getDistribution() );

    const DenseVector<ValueType>* denseX = dynamic_cast<const DenseVector<ValueType>*>( &x );
    const DenseVector<ValueType>* denseY = dynamic_cast<const DenseVector<ValueType>*>( &y );
    DenseVector<ValueType>* denseResult = dynamic_cast<DenseVector<ValueType>*>( &result );

    LAMA_ASSERT( denseX, x << ": must be DenseVector<" << Scalar::getType<ValueType>() << ">" );

    // Note: in case of beta == 0, we might skip this test

    LAMA_ASSERT( denseY, y << ": must be DenseVector<" << Scalar::getType<ValueType>() << ">" );

    LAMA_ASSERT( denseResult, result << ": must be DenseVector<" << Scalar::getType<ValueType>() << ">" );

    matrixTimesVectorImpl( *denseResult, alpha.getValue<ValueType>(), *denseX, beta.getValue<ValueType>(), *denseY );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::matrixTimesScalar( const Matrix& other, Scalar alpha )
{
    LAMA_LOG_INFO( logger, "this  = " << alpha << " * " << other );

    // should also work fine if other == *this, will not create new data

    assign( other );

    mLocalData->scale( alpha.getValue<ValueType>() );
    mHaloData->scale( alpha.getValue<ValueType>() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseMatrix<ValueType>::maxNorm() const
{
    LAMA_REGION( "Mat.Sp.maxNorm" );

    ValueType myMax = mLocalData->maxNorm();
    ValueType myMaxHalo = mHaloData->maxNorm();

    if( myMaxHalo > myMax )
    {
        myMax = myMaxHalo;
    }

    const Communicator& comm = getDistribution().getCommunicator();

    ValueType allMax = comm.max( myMax );

    LAMA_LOG_INFO( logger, "max norm: local max = " << myMax << ", global max = " << allMax );

    return allMax;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar SparseMatrix<ValueType>::maxDiffNorm( const Matrix& other ) const
{
    // Implementation works only for same row distribution, replicated col distribution
    // and the same type

    LAMA_REGION( "Mat.Sp.maxDiffNorm" );

    if( ( getDistribution() == other.getDistribution() ) && getColDistribution().isReplicated()
            && other.getColDistribution().isReplicated() && ( getValueType() == other.getValueType() ) )
    {
        const SparseMatrix<ValueType>* typedOther = dynamic_cast<const SparseMatrix<ValueType>*>( &other );
        LAMA_ASSERT_DEBUG( typedOther, "SERIOUS: wrong dynamic cast: " << other );
        return maxDiffNormImpl( *typedOther );
    }
    else if( !getColDistribution().isReplicated() )
    {
        // @todo handle maxDiffNorm on sparse matrices with column distribution
        LAMA_THROWEXCEPTION( "maxDiffNorm not available: " << *this << " has column distribution" );
    }
    else
    {
        LAMA_UNSUPPORTED( "maxDiffNorm requires temporary of " << other );
        SparseMatrix<ValueType> typedOther( other, getDistributionPtr(), getColDistributionPtr() );
        return maxDiffNormImpl( typedOther );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType SparseMatrix<ValueType>::maxDiffNormImpl( const SparseMatrix<ValueType>& other ) const
{
    // implementation only supported for same row distributions, replicated columns

    LAMA_ASSERT_EQUAL_ERROR( getDistribution(), other.getDistribution() );
    LAMA_ASSERT_ERROR( getColDistribution().isReplicated(), *this << ": not replicated column dist" );
    LAMA_ASSERT_ERROR( other.getColDistribution().isReplicated(), other << ": not replicated column dist" );

    ValueType myMaxDiff = mLocalData->maxDiffNorm( other.getLocalStorage() );

    const Communicator& comm = getDistribution().getCommunicator();

    ValueType allMaxDiff = comm.max( myMaxDiff );

    LAMA_LOG_INFO( logger, "max diff norm: local max = " << myMaxDiff << ", global max = " << allMaxDiff );

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
    return getDistribution().getCommunicator().sum( getPartitialNumValues() );
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
    const Distribution& distributionRow = getDistribution();
    const Distribution& distributionCol = getColDistribution();
    LAMA_LOG_TRACE( logger, "this(" << i << "," << j << ")" );
    ValueType myValue = 0.0;
    const IndexType iLocal = distributionRow.global2local( i );

    if( iLocal != nIndex )
    {
        LAMA_LOG_TRACE( logger, "row " << i << " is local " << iLocal );
        IndexType jLocal = distributionCol.global2local( j );

        if( nIndex != jLocal )
        {
            LAMA_LOG_TRACE( logger, "global(" << i << "," << j << ")"
                            " is local(" << iLocal << "," << jLocal << ")" );
            myValue = mLocalData->getValue( iLocal, jLocal );
            LAMA_LOG_TRACE( logger, "found local value " << myValue );
        }
        else
        {
            jLocal = mHalo.global2halo( j );
            LAMA_LOG_TRACE( logger, "global(" << i << "," << j << ")"
                            " is halo(" << iLocal << "," << jLocal << ")" );
            myValue = mHaloData->getValue( iLocal, jLocal );
            LAMA_LOG_TRACE( logger, "found halo value " << myValue );
        }
    }

    LAMA_LOG_TRACE( logger, "myValue = " << myValue );
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
    stream << getTypeName() << "( " << "size = " << mNumRows << "x" << mNumColumns << ", local = " << *mLocalData
           << ", halo = " << *mHaloData << ", row dist = " << getDistribution() << ", col dist = "
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
    if( getDistribution() != getColDistribution() )
    {
        return false;
    }

    int localDiagProperty = 0;

    if( mLocalData->hasDiagonalProperty() )
    {
        localDiagProperty = 1;
    }

    int globalDiagProperty = getDistribution().getCommunicator().min( localDiagProperty );

    return ( globalDiagProperty == 1 );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::resetDiagonalProperty()
{
    if( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "diagonal property not possible " );
    }

    this->mLocalData->resetDiagonalProperty();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
Scalar::ScalarType SparseMatrix<ValueType>::getValueType() const
{
    return Scalar::getType<ValueType>();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
std::auto_ptr<Matrix> SparseMatrix<ValueType>::create() const
{
    LAMA_LOG_INFO( logger, "SparseMatrix<ValueType>::create" );

    shared_ptr<MatrixStorage<ValueType> > newLocalData( mLocalData->create() );

    Matrix* newSparseMatrix =

        new SparseMatrix<ValueType>( newLocalData );

    // inherit the context for local and halo storage

    newSparseMatrix->setContext( mLocalData->getContextPtr(), mHaloData->getContextPtr() );

    newSparseMatrix->setCommunicationKind( this->getCommunicationKind() );

    return std::auto_ptr<Matrix>( newSparseMatrix );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
std::auto_ptr<Matrix> SparseMatrix<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "copy of " << *this );

    // copy makes deep copies of the local + halo storage, halo
    //      and uses the same distributions

    shared_ptr<MatrixStorage<ValueType> > newLocalData( mLocalData->copy() );
    shared_ptr<MatrixStorage<ValueType> > newHaloData( mHaloData->copy() );

    Matrix* newSparseMatrix =

        new SparseMatrix<ValueType>( newLocalData, newHaloData, mHalo, getDistributionPtr(), getColDistributionPtr() );

    LAMA_LOG_INFO( logger, "copy is " << *newSparseMatrix );

    return std::auto_ptr<Matrix>( newSparseMatrix );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::setIdentity()
{
    if( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( *this << ": setIdentity only supported for same row/col distribution" );
    }

    mLocalData->setIdentity( getDistributionPtr()->getLocalSize() );
    mHaloData->clear();
    mHalo.clear();

    LAMA_LOG_INFO( logger, *this << ": is no identity" );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
size_t SparseMatrix<ValueType>::getMemoryUsage() const
{
    size_t memoryUsage = mLocalData->getMemoryUsage() + mHaloData->getMemoryUsage();

    return getDistribution().getCommunicator().sum( memoryUsage );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::writeToFile(

    const std::string& fileName,
    const File::FileType fileType /* = UNFORMATTED */,
    const File::DataType dataType /* = INTERNAL */,
    const File::IndexDataType indexDataTypeIA /* = LONG */,
    const File::IndexDataType indexDataTypeJA /* = LONG */) const
{
    if( getDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file

        const Communicator& comm = getDistribution().getCommunicator();

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
        LAMA_LOG_INFO( logger, "write distributed matrix" );

        const Communicator& comm = getDistribution().getCommunicator();

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

            local.joinHalo( *mLocalData, *mHaloData, mHalo, getColDistribution() );

            local.writeToFile( comm.getSize(), comm.getRank(), fileName, fileType, dataType, indexDataTypeIA,
                               indexDataTypeJA );
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( *this << ": write to file not supported with distributions" );
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void SparseMatrix<ValueType>::readFromFile( const std::string& fileName )
{
    LAMA_REGION( "Mat.Sp.readFromFile" );

    // Take the current default communicator
    CommunicatorPtr comm = CommunicatorFactory::get();

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

/* ------------------------------------------------------------------------- */

/* ========================================================================= */

template<>
const char* SparseMatrix<float>::typeName()
{
    return "SparseMatrix<float>";
}

template<>
const char* SparseMatrix<double>::typeName()
{
    return "SparseMatrix<double>";
}

/* ------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT SparseMatrix<float> ;
template class LAMA_DLL_IMPORTEXPORT SparseMatrix<double> ;

} //namespace lama
