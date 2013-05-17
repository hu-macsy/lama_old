/**
 * @file DenseMatrix.cpp
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
 * @brief DenseMatrix.cpp
 * @author Michael Drost
 * @date 22.02.2011
 * $Id$
 */

// hpp
#include <lama/matrix/DenseMatrix.hpp>

// others
#include <lama/matrix/CSRSparseMatrix.hpp>

#include <lama/DenseVector.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/ContextAccess.hpp>
#include <lama/NoSyncToken.hpp>
#include <lama/tracing.hpp>

#include <lama/distribution/NoDistribution.hpp>
#include <lama/distribution/CyclicDistribution.hpp>
#include <lama/distribution/Redistributor.hpp>

// boost
#include <boost/scoped_array.hpp>

namespace lama
{

/* ========================================================================= */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename T>, DenseMatrix<T>::logger, "Matrix.DenseMatrix" )

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::computeOwners()
{
    // build global vector mOwners with mOwners[i] is owner of column i

    const Distribution& colDist = getColDistribution();

    LAMA_LOG_DEBUG( logger, "computerOwners for col dist = " << colDist )

    mOwners.resize( mNumColumns );

    {
        std::vector<IndexType> requiredIndexes(mNumColumns);

        for( std::vector<IndexType>::size_type i = 0; i < requiredIndexes.size(); ++i)
        {
            requiredIndexes[i] = static_cast<IndexType>( i );
        }

        if ( LAMA_LOG_TRACE_ON(logger) )
        {
            std::string s = "requiredIndexes{ ";
            for( unsigned int i = 0; i < requiredIndexes.size(); ++i)
            {
                s += " " + boost::lexical_cast<std::string>(requiredIndexes[i]);
            }
            s += " }";
            LAMA_LOG_TRACE(logger,s);
        }
        colDist.computeOwners( requiredIndexes, mOwners );
    }
    if ( LAMA_LOG_TRACE_ON(logger) )
    {
        std::string s = "mOwners{ ";
        for( std::vector<PartitionId>::size_type i = 0; i < mOwners.size(); ++i)
        {
            s += " " + boost::lexical_cast<std::string>(mOwners[i]);
        }
        s += " }";
        LAMA_LOG_TRACE(logger,s);
    }
}

/* ========================================================================= */
/*       Methods of DenseMatrix<ValueType>                                   */
/* ========================================================================= */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix()
{
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( 0, 0 ) );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const IndexType numRows, const IndexType numColumns )

    : CRTPMatrix<DenseMatrix<ValueType>,ValueType>( numRows, numColumns )
{
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( numRows, numColumns ) );
    computeOwners();
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DistributionPtr rowDist, DistributionPtr colDist )
    :

    CRTPMatrix<DenseMatrix<ValueType>,ValueType>( rowDist, colDist )
{
    computeOwners();
    allocateData();
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const DenseMatrix<ValueType>& matrix,
    DistributionPtr rowDistribution,
    DistributionPtr colDistribution )
    :

    CRTPMatrix<DenseMatrix<ValueType>,ValueType>( rowDistribution, colDistribution )
{
    LAMA_LOG_INFO( logger, "redistribute " << matrix << " to " << *this )

    if ( matrix.getColDistribution().getNumPartitions() != 1 || matrix.getDistribution().getNumPartitions() != 1 )
    {
        LAMA_THROWEXCEPTION( "Redistribution of a distributed DenseMatrix is not implemented yet" )
    }

    const Distribution& colDist = *colDistribution.get();
    const Distribution& rowDist = *rowDistribution.get();

    computeOwners();

    const IndexType numLocalRows = rowDist.getLocalSize();
    const PartitionId numColPartitions = colDist.getNumPartitions();

    LAMA_LOG_TRACE( logger, "rowDist: " << rowDist )
    LAMA_LOG_TRACE( logger, "colDist: " << colDist )
    LAMA_LOG_TRACE( logger, "colDist.comm: " << colDist.getCommunicator() )

    boost::shared_ptr<DenseStorage<ValueType> > otherData = matrix.mData[0];

    if ( rowDist.getNumPartitions() > 1 )
    {
        // only local rows of other matrix will be stored here

        otherData.reset( new DenseStorage<ValueType>( numLocalRows, getNumColumns() ) );

        localize( *otherData, *matrix.mData[0], rowDist );

        LAMA_LOG_INFO( logger, "localized for " << rowDist << ": " << *otherData )
    }

    splitColumnData( mData, *otherData, numColPartitions, mOwners );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const Matrix& other,
    DistributionPtr rowDistribution,
    DistributionPtr colDistribution )
{
    LAMA_LOG_INFO( logger, "construct copy of " << other << " for " << *this )
    assign( other );
    redistribute( rowDistribution, colDistribution );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const _MatrixStorage& other,
    DistributionPtr rowDistribution,
    DistributionPtr colDistribution )
{
    assign( other, rowDistribution, colDistribution );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> expression )
{
    Matrix::operator=( expression );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> expression )
{
    Matrix::operator=( expression );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression<Matrix,Matrix,Times> expression )
{
    Matrix::operator=( expression );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression<Scalar,Matrix,Times> expression )
{
    Matrix::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const std::string& fileName )
{
    LAMA_LOG_INFO( logger, "DenseMatrix( (fileName = " << fileName )

    boost::shared_ptr<DenseStorage<ValueType> > denseStorage( new DenseStorage<ValueType>() );

    denseStorage->readFromFile( fileName );

    LAMA_LOG_INFO( logger, "read dense storage from file " << fileName << ": " << denseStorage )

    mData.resize( 1 );
    mData[0] = denseStorage;

    Matrix::setReplicatedMatrix( denseStorage->getNumRows(), denseStorage->getNumColumns() );

    computeOwners();
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::writeToFile(

    const std::string& fileName,
    const File::FileType fileType /* = UNFORMATTED */,
    const File::DataType dataType /* = INTERNAL */,
    const File::IndexDataType indexDataTypeIA /* = LONG */,
    const File::IndexDataType indexDataTypeJA /* = LONG */) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": writeToFile( " << fileName << ", fileType = " << fileType << ", dataType = " << dataType << " )" )

    if ( getDistribution().isReplicated() && getColDistribution().isReplicated() )
    {
        // make sure that only one processor writes to file

        const Communicator& comm = getDistribution().getCommunicator();

        if ( comm.getRank() == 0 )
        {
            mData[0]->writeToFile( fileName, fileType, dataType, indexDataTypeIA, indexDataTypeJA );
        }

        // synchronization to avoid that other processors start with
        // something that might depend on the finally written file

        comm.synchronize();
    }
    else
    {
        LAMA_THROWEXCEPTION( *this << ": write to file not supported with distributions" )
    }
}

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=( const DenseMatrix<ValueType>& other )
{
    assign( other );
    return *this;
}

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=( const Matrix& other )
{
    assign( other );
    return *this;
}

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=( const Expression<Scalar,Matrix,Times> expression )
{
    Matrix::operator=( expression );
    return *this;
}

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=( const Expression<Matrix,Matrix,Times> expression )
{
    Matrix::operator=( expression );
    return *this;
}

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=(
    const Expression<Scalar,Expression<Matrix,Matrix,Times>,Times> expression )
{
    Matrix::operator=( expression );
    return *this;
}

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=(
    const Expression<Expression<Scalar,Expression<Matrix,Matrix,Times>,Times>,Expression<Scalar,Matrix,Times>,Plus> expression )
{
    Matrix::operator=( expression );
    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const DenseMatrix<ValueType>& other )
    :

    CRTPMatrix<DenseMatrix<ValueType>,ValueType>()

{
    LAMA_LOG_INFO( logger, "copy constructor( dense matrix, same value type) : " << other )
    assign( other ); // will choose the local assignment
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Matrix& other )
{
    LAMA_LOG_INFO( logger, "copy constructor( any matrix) : " << other )
    assign( other ); // will choose the local assignment
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DistributionPtr distribution )

    : CRTPMatrix<DenseMatrix<ValueType>,ValueType>( distribution, distribution )
{
    const Distribution& dist = getDistribution();

    {
        const int n = dist.getNumPartitions();
        const int numLocalRows = dist.getLocalSize();
        computeOwners();
        boost::scoped_array<int> numCols( new int[n] );
        for ( int i = 0; i < n; ++i )
        {
            numCols[i] = 0;
        }
        for ( unsigned int i = 0; i < mOwners.size(); ++i )
        {
            ++numCols[mOwners[i]];
        }

        mData.resize( n );

        LAMA_LOG_DEBUG( logger, "mData.size() = " << mData.size() )

        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( int i = 0; i < n; ++i )
        {
            //create Storage Vector
            mData[i].reset( new DenseStorage<ValueType>( numLocalRows, numCols[i] ) );
        }
        mData[0]->setIdentity();

        LAMA_LOG_DEBUG( logger, "mData[0] : " << *mData[0] << ", with data = " << mData[0]->getData() )
    }

    LAMA_LOG_INFO( logger, *this << " constructed" )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setIdentity( DistributionPtr dist )
{
    Matrix::setDistributedMatrix( dist, dist );

    computeOwners();
    allocateData();

    // Note: data is already allocated, so we just set it

    const Communicator& comm = dist->getCommunicator();

    IndexType rank = comm.getRank();
    IndexType size = comm.getSize();

    LAMA_ASSERT_EQUAL_DEBUG( size, static_cast<IndexType>( mData.size() ) )

    for ( IndexType i = 0; i < size; i++ )
    {
        if ( i == rank )
        {
            mData[i]->setIdentity();
        }
        else
        {
            mData[i]->setZero();
        }
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
bool DenseMatrix<ValueType>::isConsistent() const
{
    int consistencyErrors = 0;

    // ToDo: this implementation should use a corresponding predicate of MatrixStorage

    const IndexType numLocalRows = getDistribution().getLocalSize();

    try
    {
        Matrix::checkSettings();

        for ( size_t i = 0; i < mData.size(); ++i )
        {
            LAMA_ASSERT_EQUAL_ERROR( numLocalRows, mData[i]->getNumRows() )
            mData[i]->check( "check for consistency" );
        }
    }
    catch ( ... )
    {
        consistencyErrors = 1;
    }

    // use communicator for global reduction to make sure that all processors return same value.

    consistencyErrors = getDistribution().getCommunicator().sum( consistencyErrors );

    return 0 == consistencyErrors;
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::invert( const Matrix& other )
{
    LAMA_ASSERT_ERROR( other.getNumRows() == other.getNumColumns(),
                       "invert not allowed for non-square matrices: " << other )

    // invert supported for replicated or cyclic(n) distributed matrices

    DistributionPtr rowDist = other.getDistributionPtr();
    DistributionPtr colDist = other.getColDistributionPtr();

    DistributionPtr tmpColDist( new NoDistribution( other.getNumColumns() ) );

    if ( rowDist->isReplicated() || ( ! hasScalaPack() ) )
    {
        assign( other );
        invertReplicated();
        return;
    }

    const CyclicDistribution* cyclicDist = dynamic_cast<const CyclicDistribution*>( rowDist.get() );

    DistributionPtr tmpRowDist;

    if ( cyclicDist )
    {
        tmpRowDist = rowDist;
    }
    else
    {
        const IndexType blockSize = 64;
        CommunicatorPtr comm = rowDist->getCommunicatorPtr();
        tmpRowDist.reset( new CyclicDistribution( other.getNumRows(), blockSize, comm ) );
    }

    assign( other );
    redistribute( tmpRowDist, tmpColDist );
    invertCyclic();
    redistribute( rowDist, colDist );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
bool DenseMatrix<ValueType>::hasScalaPack()
{
    // check the LAMAInterface if ScalaPack is available ( at least on Host )

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    typename BLASInterface::SCALAPACK<ValueType>::inverse 
        inverse = loc->getInterface().BLAS.inverse<ValueType>();

    return inverse != NULL;
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::invertReplicated()
{
    LAMA_REGION( "Mat.Dense.invertReplicated" )

    DistributionPtr rowDist = getDistributionPtr();
    DistributionPtr colDist = getColDistributionPtr();

    DistributionPtr repRowDist( new NoDistribution( getNumRows() ) );
    DistributionPtr repColDist( new NoDistribution( getNumColumns() ) );

    redistribute( repRowDist, repColDist );

    // now invert the dense matrix storage

    mData[0]->invert( *mData[0] );

    redistribute( rowDist, colDist );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::invertCyclic()
{
    LAMA_REGION( "Mat.Dense.invertCyclic" )

    const Communicator& comm = getDistribution().getCommunicator();

    const Distribution& rowDist = getDistribution();

    const CyclicDistribution* cyclicDist = dynamic_cast<const CyclicDistribution*>( &rowDist );

    LAMA_ASSERT_ERROR( cyclicDist, "no cyclic distribution: " << rowDist )

    const int nb = cyclicDist->chunkSize(); // blocking factor

    ContextPtr loc = getContextPtr();  // location where inverse computation will be done

    LAMA_INTERFACE_FN_DEFAULT_T( inverse, loc, BLAS, SCALAPACK, ValueType );

    // be careful: loc might have changed to location where 'inverse' is available

    const int n = getNumRows();

    // assert square matrix

    LAMA_ASSERT_EQUAL_ERROR( getNumColumns(), n )

    DenseStorage<ValueType>& denseStorage = getLocalStorage();

    const IndexType localSize = denseStorage.getData().size();

    LAMA_ASSERT_EQUAL_ERROR( localSize, denseStorage.getNumRows() * n )

    LAMA_LOG_INFO( logger, "local dense data = " << denseStorage << ", localSize = " << localSize )

    WriteAccess<ValueType> localValues( denseStorage.getData(), loc );

    ValueType* data = localValues.get();

    LAMA_LOG_INFO( logger, "now call inverse" )

    LAMA_CONTEXT_ACCESS( loc )

    inverse( n, nb, data, comm );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::buildCSRData(
    LAMAArray<IndexType>& rowIA,
    LAMAArray<IndexType>& rowJA,
    _LAMAArray& rowValues ) const
{
    if ( getValueType() != rowValues.getValueType() )
    {
        LAMA_THROWEXCEPTION( "rowValues does not fit dense matrix type" )
    }

    LAMA_THROWEXCEPTION( "buildCSRData not available yet: ia = " << rowIA << ", ja = " << rowJA )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setCSRData(
    const LAMAArray<IndexType>& rowIA,
    const LAMAArray<IndexType>& rowJA,
    const _LAMAArray& rowValues,
    DistributionPtr,
    DistributionPtr )
{
    setCSRDataLocal( rowIA, rowJA, rowValues );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setCSRDataLocal(
    const LAMAArray<IndexType>& rowIA,
    const LAMAArray<IndexType>& rowJA,
    const _LAMAArray& rowValues ) const
{
    // build DenseStorage from the CSR data

    mData[0]->setCSRData( rowIA.size() + 1, mNumColumns, rowJA.size(), rowIA, rowJA, rowValues );

    // ToDo: split up mData[0] according to column distribution
}

template<typename ValueType>
void DenseMatrix<ValueType>::clear()
{
    Matrix::setReplicatedMatrix( 0, 0 ); // clear Matrix

    mData.resize( 1 ); // clear Data

    mData[0]->clear();
}

template<typename ValueType>
void DenseMatrix<ValueType>::allocate( const IndexType numRows, const IndexType numColumns )
{
    Matrix::setReplicatedMatrix( numRows, numColumns );

    mData.resize( 1 ); // all other storages will be freed

    LAMA_ASSERT_ERROR( mData[0], "no local data available" )

    mData[0]->allocate( mNumRows, mNumColumns );
}

template<typename ValueType>
void DenseMatrix<ValueType>::allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    LAMA_LOG_DEBUG( logger,
                    *this << " with mData[" << mData.size() << "]" << ", owners[" << mOwners.size() << "] " << " allocate row dist = " << *rowDistribution << ", col dist = " << *colDistribution )

    if ( colDistribution->isReplicated() )
    {
        mData.resize( 1 ); // all other storages will be freed

        if ( mData[0] )
        {
            // just reallocate the storage
            mData[0]->allocate( rowDistribution->getLocalSize(), colDistribution->getGlobalSize() );
        }
        else
        {
            // first time allocation
            mData[0].reset( new DenseStorage<ValueType>( rowDistribution->getLocalSize(),
                            colDistribution->getGlobalSize() ) );
        }
        if ( *colDistribution != getColDistribution() )
        {
            computeOwners();
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "coldistribution not handled here" )
    }

    Matrix::setDistributedMatrix( rowDistribution, colDistribution );
    LAMA_LOG_DEBUG( logger, *this << ": now allocated" )
}

template<typename ValueType>
void DenseMatrix<ValueType>::swap( DenseMatrix<ValueType>& other )
{
    Matrix::swapMatrix( other );

    // now swap own member variables

    std::swap( mData, other.mData );
    std::swap( mOwners, other.mOwners );
}

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const Matrix& other )
{
    LAMA_LOG_INFO( logger, "assign " << &other << " to " << this )
    LAMA_LOG_INFO( logger, "assign " << other << " to " << *this )

    if ( &other == this )
    {
        LAMA_LOG_INFO( logger, "self assign, is skpped" )
        return;
    }

    // assign will not take over sizes

    Matrix::setDistributedMatrix( other.getDistributionPtr(), other.getColDistributionPtr() );

    if ( other.getMatrixKind() == Matrix::DENSE )
    {
        LAMA_LOG_INFO( logger, "copy dense matrix" )

        switch ( other.getValueType() )
        {
        case Scalar::FLOAT:

            copyDenseMatrix( dynamic_cast<const DenseMatrix<float>&>( other ) );
            break;

        case Scalar::DOUBLE:

            copyDenseMatrix( dynamic_cast<const DenseMatrix<double>&>( other ) );
            break;

        default:
            LAMA_THROWEXCEPTION( "type of dense matrix not supported for assignment: " << other )
        }
        return;
    }

    const _SparseMatrix* sparseMatrix = dynamic_cast<const _SparseMatrix*>( &other );

    if ( sparseMatrix )
    {
        assignSparse( *sparseMatrix );
        return;
    }

    LAMA_THROWEXCEPTION( "Unsupported: assign " << other << " to " << *this )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignSparse( const _SparseMatrix& other )
{
    // @todo: this routine needs some redesign

    if ( !other.getColDistribution().isReplicated() )
    {
        DistributionPtr repColDist( new NoDistribution( other.getNumColumns() ) );

        CSRSparseMatrix<ValueType> otherCSR( other, other.getDistributionPtr(), repColDist );

        // assertion just to make sure that we do not end up in infinite recursion

        LAMA_ASSERT_DEBUG( otherCSR.getColDistribution().isReplicated(), "otherCSR not replicated columns" )

        assignSparse( otherCSR );

        splitColumns( other.getColDistributionPtr() );

        return;
    }

    // replicated columns in sparse matrix, so we can assign local data

    Matrix::setDistributedMatrix( other.getDistributionPtr(), other.getColDistributionPtr() );

    mData.resize( 1 );

    mData[0].reset( new DenseStorage<ValueType>( other.getLocalStorage() ) );

    computeOwners();
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignLocal( const _MatrixStorage& other )
{
    LAMAArray<IndexType> ia;
    LAMAArray<IndexType> ja;
    LAMAArray<ValueType> values; // get values of same type this matrix needs

    other.buildCSRData( ia, ja, values );

    setCSRDataLocal( ia, ja, values );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const _MatrixStorage& storage )
{
    LAMA_LOG_INFO( logger, "assign matrix storage = " << storage )

    const IndexType numRows = storage.getNumRows();
    const IndexType numColumns = storage.getNumColumns();

    Matrix::setReplicatedMatrix( numRows, numColumns );

    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( storage ) );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist )
{
    LAMA_LOG_INFO( logger, "assign local matrix storage = " << storage )

    Matrix::setDistributedMatrix( rowDist, colDist );

    const PartitionId numColPartitions = colDist->getNumPartitions();

    computeOwners(); // compute mapping column index -> chunk

    if ( storage.getNumRows() == rowDist->getLocalSize() )
    {
        // only format conversion of the local storage, @todo avoid it if storage is DenseStorage<ValueType>

        if ( storage.getFormat() == DENSE && storage.getValueType() == getValueType() )
        {
            const DenseStorage<ValueType>* localData = dynamic_cast<const DenseStorage<ValueType>*>( &storage );
            LAMA_ASSERT_ERROR( localData, "dynamic_cast<constDenseStorage<ValueType>*> failed: " << storage )
            splitColumnData( mData, *localData, numColPartitions, mOwners );
        }
        else if ( colDist->isReplicated() )
        {
            mData.resize( 1 );
            mData[0].reset( new DenseStorage<ValueType>( storage ) );
        }
        else
        {
            DenseStorage<ValueType> localData;
            localData.assign( storage );
            splitColumnData( mData, localData, numColPartitions, mOwners );
        }
    }
    else if ( storage.getNumRows() == rowDist->getGlobalSize() )
    {
        // we also localize the rows of the matrix

        DenseStorage<ValueType> localData;
        localData.localize( storage, *rowDist );
        splitColumnData( mData, localData, numColPartitions, mOwners );
    }
    else
    {
        LAMA_THROWEXCEPTION( storage << ": does not fit to row distribution " << *rowDist )
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::buildLocalStorage( _MatrixStorage& storage ) const
{
    if ( getColDistribution().isReplicated() )
    {
        // copy local storage with format / value conversion

        // works fine: storage.assign( *mData[0] );
        storage = *mData[0];
    }
    else
    {
        // temporary local storage with joined columns needed before

        DenseStorage<ValueType> denseStorage( getDistribution().getLocalSize(), mNumColumns );
        joinColumnData( denseStorage, mData, mOwners );
        storage = denseStorage;
    }

    LAMA_LOG_INFO( logger, "buildLocalStorage( " << *this << " ) = " << storage )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::joinColumnData(
    DenseStorage<ValueType>& result,
    const std::vector<boost::shared_ptr<DenseStorage<ValueType> > >& chunks,
    const std::vector<IndexType>& columnOwners )
{
    // Note: this is a static method, no member variables are used

    const IndexType numColumns = result.getNumColumns();
    const IndexType numRows = result.getNumRows();

    // LAMA_LOG_INFO( logger, "join column data of " << chunks.size() << " chunks to " << result )

    LAMA_ASSERT_EQUAL_ERROR( static_cast<IndexType>( columnOwners.size() ), numColumns )

    const PartitionId numColPartitions = static_cast<PartitionId>( chunks.size() );

    typedef boost::shared_ptr<HostReadAccess<ValueType> > ReadAccessPtr;

    std::vector<ReadAccessPtr> chunkRead( numColPartitions );

    // Get read access to all chunks, make some assertions for each chunk

    for ( PartitionId p = 0; p < numColPartitions; ++p )
    {
        LAMA_ASSERT_ERROR( chunks[p], "no chunk data for partition " << p )
        LAMA_ASSERT_EQUAL_ERROR( chunks[p]->getNumRows(), numRows )
        chunkRead[p].reset( new HostReadAccess<ValueType>( chunks[p]->getData() ) );
        LAMA_LOG_DEBUG( logger, "column chunk[" << p << "] : " << *chunks[p] )
    }

    std::vector<IndexType> chunkOffset( numColPartitions, 0 ); // offset for each chunk

    HostWriteAccess<ValueType> resultWrite( result.getData() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            IndexType chunkId = columnOwners[j];

            HostReadAccess<ValueType>& chunkData = *chunkRead[chunkId];

            IndexType idx = chunkOffset[chunkId]++;
            resultWrite[i * numColumns + j] = chunkData[idx];
        }
    }

    // Verify that last offset for each chunk is equal to the corresponding size

    for ( PartitionId p = 0; p < numColPartitions; ++p )
    {
        LAMA_ASSERT_EQUAL_ERROR( chunkOffset[p], chunks[p]->getNumColumns() * numRows )
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::allocateData()
{
    // mOwners are already computed, now we count them

    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, (IndexType) mOwners.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, getColDistribution().getGlobalSize() )

    const PartitionId numChunks = getColDistribution().getCommunicator().getSize();

    mData.clear();
    mData.resize( numChunks );

    const IndexType numRows = getDistribution().getLocalSize();

    if ( numChunks == 1 )
    {
        // simple case, no need to count owners for each partition

        mData[0].reset( new DenseStorage<ValueType>( numRows, mNumColumns ) );
        return;
    }

    boost::scoped_array<PartitionId> numColsPartition( new PartitionId[numChunks] );

    for ( std::vector<PartitionId>::size_type i = 0; i < mOwners.size(); ++i )
    {
        LAMA_ASSERT_DEBUG( mOwners[i] < numChunks,
                           "column owner [" << i << "] = " << mOwners[i] << " out of range, #chunks = " << numChunks )

        ++numColsPartition[mOwners[i]];
    }

    for ( PartitionId p = 0; p < numChunks; ++p )
    {
        mData[p].reset( new DenseStorage<ValueType>( numRows, numColsPartition[p] ) );
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::splitColumnData(
    std::vector<boost::shared_ptr<DenseStorage<ValueType> > >& chunks,
    const DenseStorage<ValueType>& columnData,
    const PartitionId numChunks,
    const std::vector<IndexType>& columnOwners )
{
    LAMA_LOG_INFO( logger, "split columns of " << columnData << " into " << numChunks << " chunks" )

    // Note: this is a static method, no member variables are used

    const IndexType numColumns = columnData.getNumColumns();
    const IndexType numRows = columnData.getNumRows();

    LAMA_ASSERT_EQUAL_ERROR( static_cast<IndexType>( columnOwners.size() ), numColumns )

    std::vector<PartitionId> numCols( numChunks, 0 );

    for ( std::vector<PartitionId>::size_type i = 0; i < columnOwners.size(); ++i )
    {
        LAMA_ASSERT_DEBUG( columnOwners[i] < numChunks, "owner out of range" )
        ++numCols[columnOwners[i]];
    }

    chunks.clear();
    chunks.resize( numChunks );

    typedef boost::shared_ptr<HostWriteAccess<ValueType> > WriteAccessPtr;

    std::vector<WriteAccessPtr> chunkWrite( numChunks );

    // Get write access to all chunks, make some assertions for each chunk

    for ( PartitionId p = 0; p < numChunks; ++p )
    {
        chunks[p].reset( new DenseStorage<ValueType>( numRows, numCols[p] ) );
        chunkWrite[p].reset( new HostWriteAccess<ValueType>( chunks[p]->getData() ) );
        LAMA_LOG_DEBUG( logger, "column chunk[" << p << "] : " << *chunks[p] )
    }

    std::vector<IndexType> chunkOffset( numChunks, 0 ); // offset for each chunk

    HostReadAccess<ValueType> columnDataRead( columnData.getData() );

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType j = 0; j < numColumns; ++j )
        {
            IndexType chunkId = columnOwners[j];

            HostWriteAccess<ValueType>& chunkData = *chunkWrite[chunkId];

            IndexType idx = chunkOffset[chunkId]++;
            chunkData[idx] = columnDataRead[i * numColumns + j];
        }
    }

    // Verify that last offset for each chunk is equal to the corresponding size

    for ( PartitionId p = 0; p < numChunks; ++p )
    {
        LAMA_ASSERT_EQUAL_ERROR( chunkOffset[p], numCols[p] * numRows )
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    LAMA_REGION( "Mat.Dense.redistribute" )
    if ( *rowDistribution == getDistribution() && *colDistribution == getColDistribution() )
    {
        LAMA_LOG_INFO( logger, "row and column distribtion remains unchanged" )
        return;
    }
    // Currently we only support redistribution of rows, col distribution must be replicated

    if ( getColDistribution().getNumPartitions() != 1 )
    {
        // Join all column data

        const IndexType numCols = getNumColumns();
        const IndexType numLocalRows = getDistribution().getLocalSize();

        boost::shared_ptr<DenseStorage<ValueType> > colData;
        colData.reset( new DenseStorage<ValueType>( numLocalRows, numCols ) );
        joinColumnData( *colData, mData, mOwners );

        mData.clear();
        mData.resize( 1 );
        mData[0] = colData;
        this->mColDistribution.reset( new NoDistribution( getNumColumns() ) );
    }

    redistributeRows( rowDistribution );
    splitColumns( colDistribution );
}

template<typename ValueType>
void DenseMatrix<ValueType>::splitColumns( DistributionPtr colDistribution )
{
    LAMA_ASSERT_EQUAL_ERROR( 1, getColDistribution().getNumPartitions() )

    boost::shared_ptr<DenseStorage<ValueType> > oldStorage = mData[0];

    Matrix::setDistributedMatrix( getDistributionPtr(), colDistribution );

    computeOwners(); // compute mapping column index -> chunk

    LAMA_ASSERT_EQUAL_ERROR( getDistribution().getLocalSize(), oldStorage->getNumRows() )
    LAMA_ASSERT_EQUAL_ERROR( mNumColumns, oldStorage->getNumColumns() )

    const PartitionId numColPartitions = colDistribution->getNumPartitions();

    splitColumnData( mData, *oldStorage, numColPartitions, mOwners );

    // old storage will be freed here at end of scope
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::localize(
    DenseStorage<ValueType>& local,
    const DenseStorage<ValueType>& global,
    const Distribution& rowDistribution )
{
    const IndexType numLocalRows = rowDistribution.getLocalSize();
    const IndexType numColumns = global.getNumColumns();

    LAMA_ASSERT_EQUAL_ERROR( global.getNumRows(), rowDistribution.getGlobalSize() )

    local.allocate( numLocalRows, numColumns );

    HostReadAccess<ValueType> repData( global.getData() );
    HostWriteAccess<ValueType> distData( local.getData() );

    for ( IndexType irow = 0; irow < numLocalRows; ++irow )
    {
        const IndexType globalRow = rowDistribution.local2global( irow );

        LAMA_LOG_TRACE( logger, "set local row " << irow << " with global row " << globalRow )

        for ( IndexType j = 0; j < numColumns; ++j )
        {
            distData[irow * numColumns + j] = repData[globalRow * numColumns + j];
        }
    }
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
static void replicate(
    DenseStorage<ValueType>& replicatedData,
    DenseStorage<ValueType>& distributedData,
    const Distribution& distribution )
{
    const IndexType numCols = replicatedData.getNumColumns();

    LAMA_ASSERT_EQUAL_DEBUG( numCols, distributedData.getNumColumns() )
    LAMA_ASSERT_EQUAL_DEBUG( replicatedData.getNumRows(), distribution.getGlobalSize() )
    LAMA_ASSERT_EQUAL_DEBUG( distributedData.getNumRows(), distribution.getLocalSize() )

    HostWriteAccess<ValueType> globalVals( replicatedData.getData() );
    HostReadAccess<ValueType> localVals( distributedData.getData() );

    // replicate distributed rows, each row has numCols entries

    distribution.replicateN( globalVals.get(), localVals.get(), numCols );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::redistributeRows( DistributionPtr rowDistribution )
{
    IndexType nCols = getNumColumns(); //  only global column size used here

    if ( *rowDistribution == getDistribution() )
    {
        LAMA_LOG_INFO( logger, "row distribtion remains unchanged" )
        return;
    }

    if ( rowDistribution->getNumPartitions() == 1 && getDistribution().getNumPartitions() == 1 )
    {
        LAMA_LOG_INFO( logger, "replace row distribtion, all on one processor" )
        this->setDistributionPtr( rowDistribution );
        return;
    }

    if ( getDistribution().getNumPartitions() == 1 )
    {
        DenseStorage<ValueType>& oldLocalData = *mData[0];

        // current dense matrix is replicated, we have only to assign the local part

        const IndexType numLocalRows = rowDistribution->getLocalSize();

        LAMA_LOG_INFO( logger,
                       "distribute replicated rows: use " << numLocalRows << " local rows of " << getNumRows() << " global rows" )

        DenseStorage<ValueType> newLocalData( numLocalRows, nCols );

        localize( newLocalData, oldLocalData, *rowDistribution );

        oldLocalData.swap( newLocalData );
        this->setDistributionPtr( rowDistribution );
        return;
    }

    if ( rowDistribution->getNumPartitions() == 1 )
    {
        // replicate the distributed matrix

        DenseStorage<ValueType> newLocalData( mNumRows, nCols );

        DenseStorage<ValueType>& oldLocalData = getLocalStorage();

        // replicate all rows according to the current row distribution

        replicate( newLocalData, oldLocalData, getDistribution() );

        oldLocalData.swap( newLocalData );

        this->setDistributionPtr( rowDistribution );

        return;
    }

    // So we have to reorganize data, build a Redistributor

    DenseStorage<ValueType>& oldLocalData = getLocalStorage();

    LAMA_ASSERT_EQUAL_DEBUG( nCols, oldLocalData.getNumColumns() )

    DenseStorage<ValueType> newLocalData( rowDistribution->getLocalSize(), nCols );

    Redistributor redistributor( rowDistribution, getDistributionPtr() ); // target, source distributions

    redistributor.redistributeN( newLocalData.getData(), oldLocalData.getData(), nCols );

    // LAMA_THROWEXCEPTION( "redistribution of dense rows not yet available" )

    oldLocalData.swap( newLocalData );

    this->setDistributionPtr( rowDistribution );
}

template<typename ValueType>
DenseMatrix<ValueType>::~DenseMatrix()
{
    // Note: all member variables are freed by their own destructors
}

template<typename ValueType>
void DenseMatrix<ValueType>::setContext( const ContextPtr context )
{
    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->setContext( context );
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::getRow( DenseVector<ValueType>& row, const IndexType globalRowIndex ) const
{
    LAMA_LOG_DEBUG( logger, "get values of dense matrix from row " << globalRowIndex )

    if ( !row.getDistribution().isReplicated() )
    {
        LAMA_THROWEXCEPTION( "vector for row data must be replicated" )
    }

    // on a replicated matrix each processor can fill the row

    if ( getDistribution().isReplicated() )
    {
        getLocalStorage().getRow( row.getLocalValues(), globalRowIndex );
        return;
    }

    // on a distributed matrix, owner fills row and broadcasts it

    const Communicator& comm = getDistribution().getCommunicator();

    // owner fills the row

    IndexType localRowIndex = getDistribution().global2local( globalRowIndex );

    IndexType owner = 0;

    if ( localRowIndex != nIndex )
    {
        getLocalStorage().getRow( row.getLocalValues(), localRowIndex );
        owner = comm.getRank() + 1;
        LAMA_LOG_DEBUG( logger,
                        "owner of row " << globalRowIndex << " is " << owner << ", local index = " << localRowIndex )
    }

    owner = comm.sum( owner ) - 1; // get owner via a reduction

    LAMA_ASSERT_ERROR( owner >= 0, "could not find owner of row " << globalRowIndex )

    {
        HostWriteAccess<ValueType> rowAccess( row.getLocalValues() );
        comm.bcast( rowAccess.get(), getNumColumns(), owner ); // bcast the row
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::getRow( Vector& row, const IndexType globalRowIndex ) const
{
    if ( getValueType() == row.getValueType() )
    {
        // row must be a DenseVector of same type

        DenseVector<ValueType>* typedRow = dynamic_cast<DenseVector<ValueType>*>( &row );
        LAMA_ASSERT_DEBUG( typedRow, "row is not DenseVector<Matrix::ValueType>" )
        getRow( *typedRow, globalRowIndex );
    }
    else
    {
        LAMA_THROWEXCEPTION( "Value type for row = " << row.getValueType () << " invalid"
                             ", must match type of matrix = " << getValueType() )
    }
}

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::getDiagonalImpl( DenseVector<OtherValueType>& diagonal ) const
{
    diagonal.allocate( getDistributionPtr() );

    // const cast for local storage here is safe, otherwise we have to swap

    LAMAArray<OtherValueType>& localValues = diagonal.getLocalValues();

    getLocalStorage().getDiagonal( localValues );
}

template<typename ValueType>
void DenseMatrix<ValueType>::getDiagonal( Vector& diagonal ) const
{
    if ( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    // todo: if ( diagonal.getVectorKind() == Vector::DENSE        )

    if ( true )
    {
        // Dense vector with this row distribution, so we do not need a temporary array

        if ( diagonal.getValueType() == Scalar::DOUBLE )
        {
            DenseVector<double>& denseDiagonal = dynamic_cast<DenseVector<double>&>( diagonal );
            getDiagonalImpl( denseDiagonal );
            return;
        }
        else if ( diagonal.getValueType() == Scalar::FLOAT )
        {
            DenseVector<float>& denseDiagonal = dynamic_cast<DenseVector<float>&>( diagonal );
            getDiagonalImpl( denseDiagonal );
            return;
        }
    }

    // Fallback solution with temporary arrays

    LAMAArray<ValueType> localDiagonal;
    getLocalStorage().getDiagonal( localDiagonal );
    diagonal.assign( localDiagonal, getDistributionPtr() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const Vector& diagonal )
{
    if ( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( getValueType() == Scalar::DOUBLE )
    {
        if ( typeid( diagonal ) == typeid(const DenseVector<double>&) )
        {
            const DenseVector<double>& denseDiagonal = dynamic_cast<const DenseVector<double>&>( diagonal );
            getLocalStorage().setDiagonal( denseDiagonal.getLocalValues() );
        }
        else
        {
            LAMA_THROWEXCEPTION( "Unequal data type of diagonal vector and matrix storage." )
        }
    }
    else if ( getValueType() == Scalar::FLOAT )
    {
        if ( typeid( diagonal ) == typeid(const DenseVector<float>&) )
        {
            const DenseVector<float>& denseDiagonal = dynamic_cast<const DenseVector<float>&>( diagonal );
            getLocalStorage().setDiagonal( denseDiagonal.getLocalValues() );
        }
        else
        {
            LAMA_THROWEXCEPTION( "Unequal data type of diagonal vector and matrix storage." )
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "Invalid value type of diagonal." )
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const Scalar diagonalValue )
{
    if ( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    getLocalStorage().setDiagonal( diagonalValue );
}

template<typename ValueType>
void DenseMatrix<ValueType>::scale( const Vector& vector )
{
    if ( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( getValueType() == Scalar::DOUBLE )
    {
        if ( typeid( vector ) == typeid(const DenseVector<double>&) )
        {
            const DenseVector<double>& denseVector = dynamic_cast<const DenseVector<double>&>( vector );
            getLocalStorage().scale( denseVector.getLocalValues() );
        }
        else
        {
            LAMA_THROWEXCEPTION( "Unequal data type of diagonal vector and matrix storage." )
        }
    }
    else if ( getValueType() == Scalar::FLOAT )
    {
        if ( typeid( vector ) == typeid(const DenseVector<float>&) )
        {
            const DenseVector<float>& denseVector = dynamic_cast<const DenseVector<float>&>( vector );
            getLocalStorage().scale( denseVector.getLocalValues() );
        }
        else
        {
            LAMA_THROWEXCEPTION( "Unequal data type of diagonal vector and matrix storage." )
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "Invalid value type of diagonal." )
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::scale( const Scalar scaleValue )
{
    if ( getDistribution() != getColDistribution() )
    {
        LAMA_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    getLocalStorage().scale( scaleValue );
}

template<typename ValueType>
std::vector<typename DenseMatrix<ValueType>::DenseStoragePtr>& DenseMatrix<ValueType>::getCyclicLocalValues()
{
    return mData;
}

template<typename ValueType>
const std::vector<typename DenseMatrix<ValueType>::DenseStoragePtr>& DenseMatrix<ValueType>::getCyclicLocalValues() const
{
    return mData;
}

template<typename ValueType>
Scalar DenseMatrix<ValueType>::getValue( lama::IndexType i, lama::IndexType j ) const
{
    ValueType myValue = 0.0;

    const Distribution& colDist = getColDistribution();
    const Distribution& rowDist = getDistribution();

    const Communicator& comm = rowDist.getCommunicator();

    if ( getDistribution().isLocal( i ) )
    {
        const IndexType iLocal = getDistribution().global2local( i );

        PartitionId owner = comm.getRank();

        if ( colDist.getNumPartitions() == 1 )
        {
            owner = 0;
        }

        IndexType jLocal = -1;
        if ( colDist.isLocal( j ) )
        {
            jLocal = colDist.global2local( j );
        }
        else
        {
            owner = mOwners[j];
            for ( PartitionId k = 0; k <= j; ++k )
            {
                if ( owner == mOwners[k] )
                {
                    ++jLocal;
                }
            }
        }

        LAMA_ASSERT_ERROR( jLocal != nIndex, "non local column index" )
        LAMA_LOG_TRACE( logger,
                        "getting value for index(" << i << "," << j << ")" << " which is localy ( " << iLocal << "," << jLocal << " )" )
        myValue = mData[owner]->getValue( iLocal, jLocal );
    }

    LAMA_LOG_TRACE( logger, "My value is " << myValue << " starting sum reduction to produce final result." )

    return comm.sum( myValue );
}

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesScalar( const Matrix& other, Scalar alpha )
{
    LAMA_LOG_INFO( logger, " this = " << alpha << " * " << other )

    assign( other );

    LAMA_LOG_INFO( logger, " this = other = " << *this )

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->scale( alpha.getValue<ValueType>() );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesVectorImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alphaValue,
    const DenseVector<ValueType>& denseX,
    const ValueType betaValue,
    const DenseVector<ValueType>& denseY ) const
{
    LAMA_REGION( "Mat.Dense.timesVector" )

    const LAMAArray<ValueType>& localY = denseY.getLocalValues();

    LAMAArray<ValueType>& localResult = denseResult.getLocalValues();

    ContextPtr localContext = mData[0]->getContextPtr();
    const Distribution& colDist = getColDistribution();
    const Communicator& comm = colDist.getCommunicator();
    int rank = comm.getRank();
    int n = colDist.getNumPartitions();

    mData[0]->prefetch();

    //It makes no sense to prefetch denseX because, if a transfer is started
    //the halo update needs to wait for this transfer to finish

    if ( betaValue != zero )
    {
        denseY.prefetch( localContext );
    }

    const LAMAArray<ValueType>& localX = denseX.getLocalValues();

    LAMA_LOG_INFO( logger, comm << ": matrixTimesVector, localX = " << localX << ", localY = " << localY )

    if ( n == 1 )
    {
        // replicated column distribution, only on local block, X is replicated
        // localResult = alpha * mData[0] * X + beta * localY

        const DenseStorage<ValueType>& dense = *mData[0];
        LAMA_LOG_INFO( logger, comm << ": matrixTimesVector, singe dense block = " << dense )
        dense.matrixTimesVector( localResult, alphaValue, localX, betaValue, localY );
        return;
    }

    LAMA_LOG_INFO( logger, comm << ": start pipelined multiplication." )

    int size = comm.max( localX.size() ); // largest local part of X

    mSendValues.clear();
    mReceiveValues.clear();

    LAMAArray<ValueType>* sendValues = &mSendValues;
    LAMAArray<ValueType>* recvValues = &mReceiveValues;

    const ValueType one = 1.0;

    {
        // resize the receive buffer to be big enough for largest part of X

        HostWriteOnlyAccess<ValueType> wRecvValues( *recvValues, size );

        HostWriteOnlyAccess<ValueType> wSendValues( *sendValues, size );
        HostReadAccess<ValueType> rLocalX( localX );

        // fill send buffer with local X of this processor

        int i = 0;
        for ( ; i < localX.size(); ++i )
        {
            wSendValues[i] = rLocalX[i];
        }
        for ( ; i < size; ++i )
        {
            wSendValues[i] = 0.0;
        }
    }

    if ( Matrix::ASYNCHRONOUS == Matrix::getCommunicationKind() )
    {
        LAMA_LOG_INFO( logger, comm << ": asynchronous communication" )

        // asynchronous communication always requires same sizes of arrays

        std::auto_ptr<SyncToken> st( comm.shiftAsync( *recvValues, *sendValues, 1 ) );

        LAMA_LOG_INFO( logger,
                       comm << ": matrixTimesVector, my dense block = " << *mData[rank] << ", localX = " << localX << ", localY = " << localY << ", localResult = " << localResult )

        mData[rank]->matrixTimesVector( localResult, alphaValue, localX, betaValue, localY );

        st->wait();

        // Problem: asynchronsous shift does not set correctly the size of the array recvValues

        std::swap( sendValues, recvValues );

        for ( PartitionId p = 1; p < n; ++p )
        {
            PartitionId actualPartition = comm.getNeighbor( -p );

            //handle return value to allow async communication

            if ( p < ( n - 1 ) )
            {
                st.reset( comm.shiftAsync( *recvValues, *sendValues, 1 ) );
            }
            else
            {
                st.reset( new NoSyncToken() );
            }
            LAMA_LOG_INFO( logger,
                           comm << ": matrixTimesVector, actual dense block [" << actualPartition << "] = " << *mData[actualPartition] << ", sendX = " << localX << ", localResult = " << localResult )

            // adapt the size of recvValues, that is now sendValues after swap

            LAMAArrayConstView<ValueType> x( *sendValues, 0, mData[actualPartition]->getNumColumns() );

            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, x, one, localResult );
            st->wait();
            std::swap( sendValues, recvValues );
        }
    }
    else
    {
        // for synchronous communication we can use the real needed sizes

        {
            HostWriteAccess<ValueType> wSendValues( *sendValues );
            wSendValues.resize( localX.size() );
        }

        LAMA_LOG_INFO( logger, comm << ": synchronous communication" )

        comm.shiftArray( *recvValues, *sendValues, 1 );

        // For the synchronous shift we have no problems regarding the correct sizes

        LAMA_LOG_DEBUG( logger, comm << ": send " << *sendValues << ", recv " << *recvValues )

        mData[rank]->matrixTimesVector( localResult, alphaValue, *sendValues, betaValue, localResult );

        std::swap( sendValues, recvValues );

        for ( PartitionId p = 1; p < n; ++p )
        {
            PartitionId actualPartition = comm.getNeighbor( -p );
            comm.shiftArray( *recvValues, *sendValues, 1 );
            LAMA_LOG_DEBUG( logger,
                            comm << ": send " << *sendValues << ", recv " << *recvValues << ", actual = " << actualPartition )
            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, *sendValues, one, localResult );
            std::swap( sendValues, recvValues );
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixPlusMatrix(
    const Scalar alpha,
    const Matrix& matA,
    const Scalar beta,
    const Matrix& matB )
{
    LAMA_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", A = " << matA << ", B = " << matB )

    const DenseMatrix<ValueType>* denseA = dynamic_cast<const DenseMatrix<ValueType>*>( &matA );

    LAMA_ASSERT_ERROR( denseA, "Must be dense matrix<" << getValueType() << "> : " << matA )

    const DenseMatrix<ValueType>* denseB = dynamic_cast<const DenseMatrix<ValueType>*>( &matB );

    LAMA_ASSERT_ERROR( denseB, "Must be dense matrix<" << getValueType() << "> : " << matB )

    // Now we can add sparse matrices

    matrixPlusMatrixImpl( alpha.getValue<ValueType>(), *denseA, beta.getValue<ValueType>(), *denseB );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixPlusMatrixImpl(
    const ValueType alpha,
    const DenseMatrix<ValueType>& A,
    const ValueType beta,
    const DenseMatrix<ValueType>& B )
{
    LAMA_REGION( "Mat.plusMatrix" )

    // already verified

    LAMA_ASSERT_EQUAL_DEBUG( A.getDistribution(), B.getDistribution() )
    LAMA_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getColDistribution() )

    // Now we can do it completly locally

    Matrix::setDistributedMatrix( A.getDistributionPtr(), A.getColDistributionPtr() );

    // Add matrices of each chunk

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->matrixPlusMatrix( alpha, *A.mData[i], beta, *B.mData[i] );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesMatrix(
    Matrix& result,
    const Scalar alpha,
    const Matrix& B,
    const Scalar beta,
    const Matrix& C ) const
{
    LAMA_ASSERT_ERROR( getDistribution().isReplicated(), "this->rows are distributed" )
    LAMA_ASSERT_ERROR( getColDistribution().isReplicated(), "this->cols are distributed" )
    LAMA_ASSERT_ERROR( B.getDistribution().isReplicated(), "B.rows are distributed" )
    LAMA_ASSERT_ERROR( B.getColDistribution().isReplicated(), "B.cols are distributed" )
    LAMA_ASSERT_ERROR( C.getDistribution().isReplicated(), "C.rows are distributed" )
    LAMA_ASSERT_ERROR( C.getColDistribution().isReplicated(), "C.cols are distributed" )

    // Prefetch values to the ComputeLocation

    DenseMatrix* res = dynamic_cast<DenseMatrix*>( &result );

    if ( res == NULL )
    {
        LAMA_THROWEXCEPTION( "Only DenseMatrix DenseMatrix Multiplication is supported." )
    }
    const DenseMatrix* Bp = dynamic_cast<const DenseMatrix*>( &B );
    if ( Bp == NULL )
    {
        LAMA_THROWEXCEPTION( "Only DenseMatrix DenseMatrix Multiplication is supported." )
    }
    const DenseMatrix* Cp = dynamic_cast<const DenseMatrix*>( &C );
    if ( Cp == NULL )
    {
        LAMA_THROWEXCEPTION( "Only DenseMatrix DenseMatrix Multiplication is supported." )
    }

    if ( res == this )
    {
        LAMA_LOG_DEBUG( logger, "result is aliased with this A matrix" )
    }
    else if ( res == Bp )
    {
        LAMA_LOG_DEBUG( logger, "result is aliased with B matrix" )
    }
    else if ( res == Cp && beta != 0.0 )
    {
        LAMA_LOG_DEBUG( logger, "result is aliased with C matrix" )
    }
    else
    {
        LAMA_LOG_DEBUG( logger, "result is not aliased, so allocate it correctly" )
        res->allocate( getDistributionPtr(), B.getColDistributionPtr() );
    }

    ContextPtr localContext = mData[0]->getContextPtr();
    res->prefetch( localContext );
    mData[0]->prefetch();
    Bp->prefetch( localContext );
    Cp->prefetch( localContext );

    //We are calculating with a replicated Matrix. So there is no need for an asyncronous call,
    //because we have to sync in this method anyway (returning void not SyncToken)

    // Note: any alias will be resolved by the matrix storage routine and not here
    //       as it might introduce a temporary in any case

    res->mData[0]->matrixTimesMatrix( alpha.getValue<ValueType>(), *mData[0], *Bp->mData[0], beta.getValue<ValueType>(),
                                      *Cp->mData[0] );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseMatrix<ValueType>::maxNorm() const
{
    ValueType myMaxDiff = 0.0;

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        ValueType maxDiff = mData[i]->maxNorm();

        if ( maxDiff > myMaxDiff )
        {
            myMaxDiff = maxDiff;
        }
    }

    const Communicator& comm = getDistribution().getCommunicator();

    return comm.max( myMaxDiff );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseMatrix<ValueType>::maxDiffNorm( const Matrix& other ) const
{
    // Implemenation works only for same distributions and same type

    if ( ( getDistribution() == other.getDistribution() ) && ( getColDistribution() == other.getColDistribution() )
            && ( getValueType() == other.getValueType() ) )
    {
        const DenseMatrix<ValueType>* typedOther = dynamic_cast<const DenseMatrix<ValueType>*>( &other );
        LAMA_ASSERT_DEBUG( typedOther, "SERIOUS: wrong dynamic cast: " << other )
        return maxDiffNormImpl( *typedOther );
    }
    else
    {
        LAMA_UNSUPPORTED( "maxDiffNorm requires temporary of " << other )
        DenseMatrix<ValueType> typedOther( other, getDistributionPtr(), getColDistributionPtr() );
        return maxDiffNormImpl( typedOther );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseMatrix<ValueType>::maxDiffNormImpl( const DenseMatrix<ValueType>& other ) const
{
    // implementation only supported for same distributions

    LAMA_ASSERT_EQUAL_ERROR( getDistribution(), other.getDistribution() )
    LAMA_ASSERT_EQUAL_ERROR( getColDistribution(), other.getColDistribution() )

    ValueType myMaxDiff = 0.0;

    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        ValueType maxDiff = mData[i]->maxDiffNorm( *other.mData[i] );
        if ( maxDiff > myMaxDiff )
        {
            myMaxDiff = maxDiff;
        }
    }

    const Communicator& comm = getDistribution().getCommunicator();

    return comm.max( myMaxDiff );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::prefetch() const
{
    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        mData[i]->prefetch();
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::prefetch( lama::ContextPtr loc ) const
{
    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        mData[i]->prefetch( loc );
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::wait() const
{
    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        mData[i]->wait();
    }
}

template<typename ValueType>
const DenseStorage<ValueType>& DenseMatrix<ValueType>::getLocalStorage() const
{
    LAMA_ASSERT_ERROR( mData.size() > 0, "no local values allocated" )

    if ( mData.size() == 1 )
    {
        return *mData[0];
    }

    LAMA_ASSERT_EQUAL_ERROR( getDistribution(), getColDistribution() )

    const PartitionId myRank = getDistribution().getCommunicator().getRank();

    return *mData[myRank];
}

template<typename ValueType>
DenseStorage<ValueType>& DenseMatrix<ValueType>::getLocalStorage()
{
    LAMA_ASSERT_ERROR( mData.size() > 0, "no local values allocated" )

    if ( mData.size() == 1 )
    {
        return *mData[0];
    }

    LAMA_ASSERT_EQUAL_ERROR( getDistribution(), getColDistribution() )

    const PartitionId myRank = getDistribution().getCommunicator().getRank();

    return *mData[myRank];
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getLocalNumValues() const
{
    // only locally stored number of values
    return getDistribution().getLocalSize() * mNumColumns;
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getLocalNumRows() const
{
    // only locally stored number of values
    return getDistribution().getLocalSize();
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getLocalNumColumns() const
{
    // only locally stored number of values
    return getColDistribution().getLocalSize();
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getNumValues() const
{
    IndexType myNumValues = 0;

    for ( size_t k = 0; k < mData.size(); ++k )
    {
        myNumValues += mData[k]->getNumValues();
    }

    return getDistribution().getCommunicator().sum( myNumValues );
}

template<typename ValueType>
bool DenseMatrix<ValueType>::hasDiagonalProperty() const
{
    // just a dummy
    return false;
}

template<typename ValueType>
void DenseMatrix<ValueType>::resetDiagonalProperty()
{
    // just a dummy
}

template<typename ValueType>
void DenseMatrix<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DenseMatrix<" << Scalar::getType<ValueType>() << ">(" << mNumRows << "x" << mNumColumns << ", rowdist = "
           << getDistribution() << ", coldist = " << getColDistribution() << ")";
}

template<typename ValueType>
Scalar::ScalarType DenseMatrix<ValueType>::getValueType() const
{
    return Scalar::getType<ValueType>();
}

template<typename ValueType>
DenseMatrix<ValueType>* DenseMatrix<ValueType>::create() const
{
    LAMA_LOG_INFO( logger, "DenseMatrix<ValueType>::create" )

    return new DenseMatrix<ValueType>();
}

template<typename ValueType>
DenseMatrix<ValueType>* DenseMatrix<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "DenseMatrix<ValueType>::copy" )

    return  new DenseMatrix<ValueType>( *this );
}

template<typename ValueType>
size_t DenseMatrix<ValueType>::getMemoryUsage() const
{
    size_t memoryUsage = 0;

    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        memoryUsage += mData[i]->getMemoryUsage();
    }

    return getDistribution().getCommunicator().sum( memoryUsage );
}

/* ========================================================================= */

template<typename ValueType>
const char* DenseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

/* ========================================================================= */

template<>
const char* DenseMatrix<float>::typeName()
{
    return "DenseMatrix<float>";
}

template<>
const char* DenseMatrix<double>::typeName()
{
    return "DenseMatrix<double>";
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

template class LAMA_DLL_IMPORTEXPORT DenseMatrix<float> ;
template class LAMA_DLL_IMPORTEXPORT DenseMatrix<double> ;

}
