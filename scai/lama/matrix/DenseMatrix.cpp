/**
 * @file DenseMatrix.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief DenseMatrix.cpp
 * @author Michael Drost
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/matrix/DenseMatrix.hpp>

// local library
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/mepr/DenseMatrixWrapper.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/unique_ptr.hpp>
#include <scai/common/ScalarType.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::hmemo;
using namespace scai::dmemo;

namespace scai
{

using common::unique_ptr;
using common::TypeTraits;
using common::scoped_array;
using utilskernel::LAMAKernel;

namespace lama
{

/* ========================================================================= */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseMatrix<ValueType>::logger, "Matrix.DenseMatrix" )

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::computeOwners()
{
    // build vector mOwners with mOwners[i] is owner of column i
    // Note: this vector is replicated for all processors, has globalSize entries

    const Distribution& colDist = getColDistribution();

    SCAI_LOG_DEBUG( logger, "computerOwners for col dist = " << colDist )

    // Note: colDist.globalSize() == mNumColumns 

    HArray<IndexType> indexes;   // will contain all column indexes to get all owners

    utilskernel::HArrayUtils::setOrder( indexes, mNumColumns );

    colDist.computeOwners( mOwners, indexes );

    SCAI_ASSERT_EQ_DEBUG( mNumColumns, mOwners.size(), "Serious mismatch, probably due to wrong distribution" );
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

    : CRTPMatrix<DenseMatrix<ValueType>, ValueType>( numRows, numColumns )
{
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( numRows, numColumns ) );
    computeOwners();
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DistributionPtr rowDist, DistributionPtr colDist )
    :

    CRTPMatrix<DenseMatrix<ValueType>, ValueType>( rowDist, colDist )
{
    computeOwners();
    allocateData();   // will initialize it with zero
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const DenseMatrix<ValueType>& other,
    DistributionPtr rowDistribution,
    DistributionPtr colDistribution )
{
    // just do the same as with any arbitrary matrix
    SCAI_LOG_INFO( logger, "construct copy of " << other << " for " << *this )
    assign( other );
    redistribute( rowDistribution, colDistribution );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    const Matrix& other,
    DistributionPtr rowDistribution,
    DistributionPtr colDistribution )
{
    SCAI_LOG_INFO( logger, "construct copy of " << other << " for " << *this )
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
DenseMatrix<ValueType>::DenseMatrix( const Expression_SMM_SM& expression )
{
    Matrix::operator=( expression );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SMM& expression )
{
    Matrix::operator=( expression );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SM_SM& expression )
{
    Matrix::operator=( expression );
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SM& expression )
{
    Matrix::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const std::string& fileName )
{
    SCAI_LOG_INFO( logger, "DenseMatrix( (fileName = " << fileName )

    // initialize mData, so getLocalStorage() returns valid data

    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( 0, 0 ) );

    this->readFromFile( fileName );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=( const DenseMatrix<ValueType>& other )
{
    // override the default assignment operator
    assign( other );
    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const DenseMatrix<ValueType>& other )
    :

    CRTPMatrix<DenseMatrix<ValueType>, ValueType>()

{
    SCAI_LOG_INFO( logger, "copy constructor( dense matrix, same value type) : " << other )
    assign( other ); // will choose the local assignment
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Matrix& other )
{
    SCAI_LOG_INFO( logger, "copy constructor( any matrix) : " << other )
    assign( other ); // will choose the local assignment
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DistributionPtr distribution )

    : CRTPMatrix<DenseMatrix<ValueType>, ValueType>( distribution, distribution )
{
    const Distribution& dist = getRowDistribution();
    {
        const int nPartitions = dist.getNumPartitions();

        const int numLocalRows = dist.getLocalSize();

        computeOwners();
       
        utilskernel::LArray<IndexType> numCols;

        utilskernel::HArrayUtils::bucketCount( numCols, mOwners, nPartitions );

        SCAI_ASSERT_EQUAL( numCols.sum(), mOwners.size(), "serious mismatch due to illegal owner" )

        mData.resize( nPartitions );  

        ReadAccess<IndexType> rNumCols( numCols );

        for ( int i = 0; i < nPartitions; ++i )
        {
            //create Storage Vector

            mData[i].reset( new DenseStorage<ValueType>( numLocalRows, rNumCols[i] ) );
        }

        mData[0]->setDiagonalImpl( ValueType( 1 ) );

        SCAI_LOG_DEBUG( logger, "mData[0] : " << *mData[0] << ", with data = " << mData[0]->getData() )
    }

    SCAI_LOG_INFO( logger, *this << " constructed" )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setIdentity( DistributionPtr dist )
{
    Matrix::setDistributedMatrix( dist, dist );
    computeOwners();
    allocateData();   // initialized with zero
    // Note: data is already allocated, so we just set it
    const Communicator& comm = dist->getCommunicator();
    IndexType rank = comm.getRank();
    IndexType size = comm.getSize();
    SCAI_ASSERT_EQUAL_DEBUG( size, static_cast<IndexType>( mData.size() ) )

    for ( IndexType i = 0; i < size; i++ )
    {
        SCAI_LOG_INFO( logger, "identity, mData[" << i << "] = " << *mData[i] );

        if ( i == rank )
        {
            mData[i]->setDiagonal( ValueType( 1 ) );
        }
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setDenseData(
    DistributionPtr rowDist,
    DistributionPtr colDist,
    const _HArray& values,
    const Scalar eps )
{
    DistributionPtr tmpReplicatedColDistribution = colDist;
    const IndexType n = rowDist->getLocalSize();
    const IndexType m = colDist->getGlobalSize();

    // splitting of the column data will be done after setting full column data

    if ( !colDist->isReplicated() )
    {
        tmpReplicatedColDistribution.reset( new NoDistribution( m ) );
    }

    Matrix::setDistributedMatrix( rowDist, tmpReplicatedColDistribution );
    // due to temporary replicated col distribution, mData has only one entry
    mData[0]->setDenseData( n, m, values, eps.getValue<ValueType>() );
    SCAI_LOG_INFO( logger,
                   "Dense matrix, row dist = " << *rowDist << " filled locally with " << ( n * m ) << " values, now split for col dist = " << *colDist );

    if ( !colDist->isReplicated() )
    {
        splitColumns( colDist );

        for ( int i = 0; i < colDist->getCommunicator().getSize(); ++i )
        {
            SCAI_LOG_DEBUG( logger, "mData[" << i << "] = " << *mData[i] );
        }
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setCSRData(
    DistributionPtr rowDist,
    DistributionPtr colDist,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    DistributionPtr tmpReplicatedColDistribution = colDist;
    const IndexType n = rowDist->getLocalSize();
    const IndexType m = colDist->getGlobalSize();

    // splitting of the column data will be done after setting full column data

    if ( !colDist->isReplicated() )
    {
        tmpReplicatedColDistribution.reset( new NoDistribution( m ) );
    }

    Matrix::setDistributedMatrix( rowDist, tmpReplicatedColDistribution );
    // due to temporary replicated col distribution, mData has only one entry
    mData[0]->setCSRData( n, m, numValues, ia, ja, values );

    if ( !colDist->isReplicated() )
    {
        splitColumns( colDist );
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setDIAData(
    DistributionPtr /*rowDist*/,
    DistributionPtr /*colDist*/,
    const IndexType /*numDiagonals*/,
    const HArray<IndexType>& /*offsets*/,
    const _HArray& /*values*/ )
{
    COMMON_THROWEXCEPTION ( "not yet implemented" )

    // from setCSRData

    // DistributionPtr tmpReplicatedColDistribution = colDist;
    // const IndexType n = rowDist->getLocalSize();
    // const IndexType m = colDist->getGlobalSize();

    // // splitting of the column data will be done after setting full column data

    // if ( !colDist->isReplicated() )
    // {
    //     tmpReplicatedColDistribution.reset( new NoDistribution( m ) );
    // }

    // Matrix::setDistributedMatrix( rowDist, tmpReplicatedColDistribution );
    // // due to temporary replicated col distribution, mData has only one entry
    // mData[0]->setCSRData( n, m, numValues, ia, ja, values );

    // if ( !colDist->isReplicated() )
    // {
    //     splitColumns( colDist );
    // }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
bool DenseMatrix<ValueType>::isConsistent() const
{
    IndexType consistencyErrors = 0;
    // ToDo: this implementation should use a corresponding predicate of MatrixStorage
    const IndexType numLocalRows = getRowDistribution().getLocalSize();

    try
    {
        Matrix::checkSettings();

        for ( size_t i = 0; i < mData.size(); ++i )
        {
            SCAI_ASSERT_EQUAL_ERROR( numLocalRows, mData[i]->getNumRows() )
            mData[i]->check( "check for consistency" );
        }
    }
    catch ( ... )
    {
        consistencyErrors = 1;
    }

    // use communicator for global reduction to make sure that all processors return same value.
    consistencyErrors = getRowDistribution().getCommunicator().sum( consistencyErrors );
    return 0 == consistencyErrors;
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::invert( const Matrix& other )
{
    SCAI_ASSERT_ERROR( other.getNumRows() == other.getNumColumns(),
                       "invert not allowed for non-square matrices: " << other )
    // invert supported for replicated or cyclic(n) distributed matrices
    DistributionPtr rowDist = other.getRowDistributionPtr();
    DistributionPtr colDist = other.getColDistributionPtr();
    DistributionPtr tmpColDist( new NoDistribution( other.getNumColumns() ) );

    assign( other );
    invertReplicated();
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::invertReplicated()
{
    SCAI_REGION( "Mat.Dense.invertReplicated" )
    DistributionPtr rowDist = getRowDistributionPtr();
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
void DenseMatrix<ValueType>::buildCSRData(
    HArray<IndexType>& rowIA,
    HArray<IndexType>& rowJA,
    _HArray& rowValues ) const
{
    if ( getValueType() != rowValues.getValueType() )
    {
        COMMON_THROWEXCEPTION( "rowValues does not fit dense matrix type" )
    }

    COMMON_THROWEXCEPTION( "buildCSRData not available yet: ia = " << rowIA << ", ja = " << rowJA )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setCSRData(
    const HArray<IndexType>& rowIA,
    const HArray<IndexType>& rowJA,
    const _HArray& rowValues,
    DistributionPtr,
    DistributionPtr )
{
    setCSRDataLocal( rowIA, rowJA, rowValues );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setCSRDataLocal(
    const HArray<IndexType>& rowIA,
    const HArray<IndexType>& rowJA,
    const _HArray& rowValues ) const
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
    SCAI_ASSERT_ERROR( mData[0], "no local data available" )
    mData[0]->allocate( mNumRows, mNumColumns );
}

template<typename ValueType>
void DenseMatrix<ValueType>::allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    SCAI_LOG_DEBUG( logger,
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
        COMMON_THROWEXCEPTION( "coldistribution not handled here" )
    }

    Matrix::setDistributedMatrix( rowDistribution, colDistribution );
    SCAI_LOG_DEBUG( logger, *this << ": now allocated" )
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
void DenseMatrix<ValueType>::assignTranspose( const Matrix& other  )
{
    SCAI_LOG_INFO( logger, "assign transposed " << other << " to " << *this )
    const DenseMatrix<ValueType>* denseMatrix = dynamic_cast<const DenseMatrix<ValueType>*>( &other );

    if ( denseMatrix )
    {
        assignTransposeImpl( *denseMatrix );
    }
    else
    {
        COMMON_THROWEXCEPTION( "DenseMatrix::assignTranspose currently only implemented for dense matrices of same type" )
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::assignTransposeImpl( const DenseMatrix<ValueType>& Mat )
{
    const Communicator& comm = Mat.getRowDistribution().getCommunicator();
    IndexType size = comm.getSize();
    DistributionPtr distRow = Mat.getRowDistributionPtr();
    DistributionPtr distCol = Mat.getColDistributionPtr();

    if ( size == 1 )        // localTranspose == globalTranpose, if processor nr == 1
    {
        if ( this != &Mat )
        {
            assign( Mat );
        }

        mData[0]->transposeImpl();
        redistribute( distCol, distRow );
    }
    else
    {
        //new storage, distribution already changed
        DenseMatrix<ValueType> targetMat( distCol, distRow );

        //local transpose of Mat
        for ( IndexType i = 0; i < size; ++i )
        {
            Mat.mData[i]->transposeImpl();
        }

        //preparation for mpi all2allv
        IndexType* receiveSizes = new IndexType[size];
        ValueType** recvBuffer = new ValueType*[size];

        for ( IndexType i = 0; i < size; ++i )
        {
            IndexType localSize = targetMat.mData[i]->getData().size();
            recvBuffer[i] = new ValueType[localSize];
            WriteAccess<ValueType> wData( targetMat.mData[i]->getData() );
            //local data
            recvBuffer[i] = wData.get();
            //local data sizes
            receiveSizes[i] = localSize;
        }

        IndexType* sendSizes = new IndexType[size];
        ValueType** sendBuffer = new ValueType*[size];

        for ( IndexType i = 0; i < size; ++i )
        {
            IndexType localSize =  Mat.mData[i]->getData().size();
            sendBuffer[i] = new ValueType[ localSize];
            WriteAccess<ValueType> wDatas( Mat.mData[i]->getData() );
            //local datal
            sendBuffer[i] = wDatas.get();
            //local data sizes
            sendSizes[i] = localSize;
        }

        //MPI call
        comm.all2allv( recvBuffer, receiveSizes, sendBuffer, sendSizes );

        //transpose back of Mat (A^t)^t = A
        if ( this != &Mat ) // no need if we override Mat anyways
        {
            for ( IndexType i = 0; i < size; ++i )
            {
                Mat.mData[i]->transposeImpl();
            }
        }

        *this = targetMat;
        delete [] sendSizes;
        delete [] receiveSizes;
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const Matrix& other )
{
    SCAI_LOG_INFO( logger, "assign " << &other << " to " << this )
    SCAI_LOG_INFO( logger, "assign " << other << " to " << *this )

    if ( &other == this )
    {
        SCAI_LOG_INFO( logger, "self assign, is skpped" )
        return;
    }

    // assign will not take over sizes
    Matrix::setDistributedMatrix( other.getRowDistributionPtr(), other.getColDistributionPtr() );

    if ( other.getMatrixKind() == Matrix::DENSE )
    {
        SCAI_LOG_INFO( logger, "copy dense matrix" )
        mepr::DenseMatrixWrapper<ValueType, SCAI_ARITHMETIC_HOST_LIST>::assignDenseImpl( *this, other );
        return;
    }
    else if ( other.getMatrixKind() == Matrix::SPARSE )
    {
        SCAI_LOG_INFO( logger, "copy sparse matrix" )
        mepr::DenseMatrixWrapper<ValueType, SCAI_ARITHMETIC_HOST_LIST>::assignSparseImpl( *this, other );
        return;
    }

    SCAI_LOG_TRACE( logger, "Unsupported assign" )
    COMMON_THROWEXCEPTION( "Unsupported: assign " << other << " to " << *this )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignSparse( const CRTPMatrix<SparseMatrix<ValueType>, ValueType>& other )
{
// @todo: this routine needs some redesign
    if ( !other.getColDistribution().isReplicated() )
    {
        DistributionPtr repColDist( new NoDistribution( other.getNumColumns() ) );
        CSRSparseMatrix<ValueType> otherCSR( other, other.getRowDistributionPtr(), repColDist );
// assertion just to make sure that we do not end up in infinite recursion
        SCAI_ASSERT_DEBUG( otherCSR.getColDistribution().isReplicated(), "otherCSR not replicated columns" )
        assignSparse( otherCSR );
        splitColumns( other.getColDistributionPtr() );
        return;
    }

// replicated columns in sparse matrix, so we can assign local data
    Matrix::setDistributedMatrix( other.getRowDistributionPtr(), other.getColDistributionPtr() );
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( other.getLocalStorage() ) );
    computeOwners();
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignLocal( const _MatrixStorage& other )
{
    HArray<IndexType> ia;
    HArray<IndexType> ja;
    HArray<ValueType> values; // get values of same type this matrix needs
    other.buildCSRData( ia, ja, values );
    setCSRDataLocal( ia, ja, values );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const _MatrixStorage& storage )
{
    SCAI_LOG_INFO( logger, "assign matrix storage = " << storage )
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
    SCAI_LOG_INFO( logger, "assign local matrix storage = " << storage )
    Matrix::setDistributedMatrix( rowDist, colDist );
    const PartitionId numColPartitions = colDist->getNumPartitions();
    computeOwners(); // compute mapping column index -> chunk

    if ( storage.getNumRows() == rowDist->getLocalSize() )
    {
        // only format conversion of the local storage, @todo avoid it if storage is DenseStorage<ValueType>

        if ( storage.getFormat() == Format::DENSE && storage.getValueType() == getValueType() )
        {
            const DenseStorage<ValueType>* localData = dynamic_cast<const DenseStorage<ValueType>*>( &storage );
            SCAI_ASSERT_ERROR( localData, "dynamic_cast<constDenseStorage<ValueType>*> failed: " << storage )

            splitColumnData( mData, *localData, numColPartitions, mOwners );
        }
        else if ( colDist->isReplicated() )
        {
            common::shared_ptr<DenseStorage<ValueType> > dataPtr( new DenseStorage<ValueType>( storage ) );

            mData.resize( 1 );
            mData[0] = dataPtr;
        }
        else
        {
            DenseStorage<ValueType> localData( storage );
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
        COMMON_THROWEXCEPTION( storage << ": does not fit to row distribution " << *rowDist )
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
        const IndexType numLocalRows = getRowDistribution().getLocalSize();
        DenseStorage<ValueType> denseStorage( numLocalRows, mNumColumns );
        joinColumnData( denseStorage.getData(), 0, numLocalRows );
        storage = denseStorage;
    }

    SCAI_LOG_INFO( logger, "buildLocalStorage( " << *this << " ) = " << storage )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::joinColumnData(
    HArray<ValueType>& result,
    const IndexType firstRow,
    const IndexType nRows ) const
{
    SCAI_LOG_DEBUG( logger, "join column data, firstRow = " << firstRow << ", nRows = " << nRows << ", result = " << result )
    const IndexType ncol = getNumColumns();
    // make sure that the array mOwners is set correctly
    SCAI_ASSERT_EQUAL_ERROR( static_cast<IndexType>( mOwners.size() ), ncol )

    const PartitionId numColPartitions = static_cast<PartitionId>( mData.size() );

    ContextPtr hostContext = Context::getHostPtr();

    // Get read access to all chunks, make some assertions for each chunk

    IndexType numGlobalColumns = 0;

    typedef const ValueType* PtrType;

    common::scoped_array<PtrType> chunkPtr( new PtrType [ numColPartitions ] );

    for ( PartitionId p = 0; p < numColPartitions; ++p )
    {
        SCAI_ASSERT_ERROR( mData[p], "no chunk data for partition " << p )

        IndexType numLocalColumns = mData[p]->getNumColumns();

        ReadAccess<ValueType> chunkRead( mData[p]->getData(), hostContext );

        chunkPtr[p] = chunkRead.get() + numLocalColumns * firstRow;

        numGlobalColumns += numLocalColumns;

        SCAI_LOG_DEBUG( logger, "column chunk[" << p << "] : " << *mData[p] )
    }

    SCAI_ASSERT_EQUAL_ERROR( numGlobalColumns, ncol );   // local column sizes must add to global column size

    SCAI_LOG_DEBUG( logger, "resize result to " << ncol << " x " << nRows )

    WriteOnlyAccess<ValueType> wResult( result, hostContext, ncol * nRows );

    ReadAccess<PartitionId> rOwners( mOwners, hostContext );

    // gather data of chunks, each chunk data is contiguous

    IndexType pos = 0;

    // No OpenMP parallelization, would require chunkPtr array for each thread

    for ( IndexType i = firstRow; i < firstRow + nRows ; ++i )
    {
        for ( IndexType j = 0; j < ncol; ++j )
        {
            IndexType chunkId = rOwners[j];
            wResult[ pos++ ] = *chunkPtr[chunkId]++;
        }
    }

    SCAI_LOG_DEBUG( logger, "ready join column data" )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::allocateData()
{
// mOwners are already computed, now we count them
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, ( IndexType ) mOwners.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, getColDistribution().getGlobalSize() )
    const PartitionId numChunks = getColDistribution().getCommunicator().getSize();
    mData.clear();
    mData.resize( numChunks );
    const IndexType numRows = getRowDistribution().getLocalSize();
    SCAI_LOG_INFO( logger, "build " << numChunks << " data arrays for numRows = " << numRows );

    if ( numChunks == 1 )
    {
        // simple case, no need to count owners for each partition

        mData[0].reset( new DenseStorage<ValueType>( numRows, mNumColumns ) );
        return;
    }

    utilskernel::LArray<IndexType> numColsPartition;
 
    utilskernel::HArrayUtils::bucketCount( numColsPartition, mOwners, numChunks );

    ContextPtr ctx = Context::getHostPtr();

    ReadAccess<IndexType> rSizes( numColsPartition, ctx );

    IndexType count = 0; // sum up the sizes, verify correct sum

    for ( PartitionId p = 0; p < numChunks; ++p )
    {
        count += rSizes[p];
        mData[p].reset( new DenseStorage<ValueType>( numRows, rSizes[p] ) );
    }

    SCAI_ASSERT_EQ_ERROR( count, mNumColumns, "Illegal owners." )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::splitColumnData(
    std::vector<common::shared_ptr<DenseStorage<ValueType> > >& chunks,
    const DenseStorage<ValueType>& columnData,
    const PartitionId numPartitions,
    const HArray<PartitionId>& columnOwners )
{
    if ( numPartitions == 1 )
    {
        SCAI_LOG_INFO( logger, "split columns of " << columnData << " unnecessary, #partitions = 1" )

        chunks.resize( 1 );

        if ( !chunks[0].get() ) 
        {
            chunks[0].reset( new DenseStorage<ValueType>( columnData ) );
        }
        else
        {
            // can be aliased here, avoids copy

            *chunks[0] = columnData;
        }

        return;
    }

    SCAI_LOG_INFO( logger, "split columns of " << columnData << " into " << numPartitions << " partition chunks" )

    const IndexType numColumns = columnData.getNumColumns();
    const IndexType numRows = columnData.getNumRows();

    SCAI_ASSERT_EQUAL_ERROR( columnOwners.size(), numColumns )

    utilskernel::LArray<IndexType> offsets;
    utilskernel::LArray<IndexType> perm;

    utilskernel::HArrayUtils::bucketSort( offsets, perm, columnOwners, numPartitions );

    // offsets array, last entry stands for number of elements sorted into buckets

    SCAI_ASSERT_EQ_DEBUG( numColumns, offsets[numPartitions], "Illegal column owners" )
    SCAI_ASSERT_EQ_DEBUG( columnOwners.size(), perm.size(), "Illegal column owners" )

    chunks.clear();
    chunks.resize( numPartitions );
  
    ContextPtr ctx = Context::getHostPtr();

    ReadAccess<ValueType> columnDataRead( columnData.getData(), ctx );
    ReadAccess<IndexType> rPerm( perm, ctx );
    ReadAccess<IndexType> rOffsets( offsets, ctx );
 
    for ( PartitionId p = 0; p < numPartitions; ++p )
    {
        IndexType lb = rOffsets[p]; 
        IndexType ub = rOffsets[p+1];

        chunks[p].reset( new DenseStorage<ValueType>( numRows, ub - lb ) );

        WriteAccess<ValueType> wChunkData( chunks[p]->getData(), ctx );
   
        IndexType pos = 0;  // traversing the elements of chunk data for p-th partition

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType j = lb; j < ub; ++j )
            {
                IndexType jj = rPerm[j]; 
                wChunkData[pos++] = columnDataRead[ i * numColumns + jj ];
            }
        }
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::redistribute( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    SCAI_REGION( "Mat.Dense.redistribute" )

    if ( *rowDistribution == getRowDistribution() && *colDistribution == getColDistribution() )
    {
        SCAI_LOG_INFO( logger, "row and column distribtion remains unchanged" )
        return;
    }

// Currently we only support redistribution of rows, col distribution must be replicated

    if ( getColDistribution().getNumPartitions() != 1 )
    {
// Join all column data
        const IndexType numCols = getNumColumns();
        const IndexType numLocalRows = getRowDistribution().getLocalSize();
        common::shared_ptr<DenseStorage<ValueType> > colData;
        colData.reset( new DenseStorage<ValueType>( numLocalRows, numCols ) );
        joinColumnData( colData->getData(), 0, numLocalRows );
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
    SCAI_ASSERT_EQUAL_ERROR( 1, getColDistribution().getNumPartitions() )
    common::shared_ptr<DenseStorage<ValueType> > oldStorage = mData[0];
    Matrix::setDistributedMatrix( getRowDistributionPtr(), colDistribution );
    computeOwners(); // compute mapping column index -> chunk
    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution().getLocalSize(), oldStorage->getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, oldStorage->getNumColumns() )
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
    SCAI_ASSERT_EQUAL_ERROR( global.getNumRows(), rowDistribution.getGlobalSize() )
    local.allocate( numLocalRows, numColumns );
    ContextPtr contextPtr = Context::getHostPtr();
    ReadAccess<ValueType> repData( global.getData(), contextPtr );
    WriteAccess<ValueType> distData( local.getData(), contextPtr );

    for ( IndexType irow = 0; irow < numLocalRows; ++irow )
    {
        const IndexType globalRow = rowDistribution.local2global( irow );
        SCAI_LOG_TRACE( logger, "set local row " << irow << " with global row " << globalRow )

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
    SCAI_ASSERT_EQUAL_DEBUG( numCols, distributedData.getNumColumns() )
    SCAI_ASSERT_EQUAL_DEBUG( replicatedData.getNumRows(), distribution.getGlobalSize() )
    SCAI_ASSERT_EQUAL_DEBUG( distributedData.getNumRows(), distribution.getLocalSize() )
    ContextPtr contextPtr = Context::getHostPtr();
    WriteAccess<ValueType> globalVals( replicatedData.getData(), contextPtr );
    ReadAccess<ValueType> localVals( distributedData.getData(), contextPtr );
// replicate distributed rows, each row has numCols entries
    distribution.replicateN( globalVals.get(), localVals.get(), numCols );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::redistributeRows( DistributionPtr rowDistribution )
{
    IndexType nCols = getNumColumns(); //  only global column size used here

    if ( *rowDistribution == getRowDistribution() )
    {
        SCAI_LOG_INFO( logger, "row distribtion remains unchanged" )
        return;
    }

    if ( rowDistribution->getNumPartitions() == 1 && getRowDistribution().getNumPartitions() == 1 )
    {
        SCAI_LOG_INFO( logger, "replace row distribtion, all on one processor" )
        this->setDistributionPtr( rowDistribution );
        return;
    }

    if ( getRowDistribution().getNumPartitions() == 1 )
    {
        DenseStorage<ValueType>& oldLocalData = *mData[0];
// current dense matrix is replicated, we have only to assign the local part
        const IndexType numLocalRows = rowDistribution->getLocalSize();
        SCAI_LOG_INFO( logger,
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
        replicate( newLocalData, oldLocalData, getRowDistribution() );
        oldLocalData.swap( newLocalData );
        this->setDistributionPtr( rowDistribution );
        return;
    }

// So we have to reorganize data, build a Redistributor
    DenseStorage<ValueType>& oldLocalData = getLocalStorage();
    SCAI_ASSERT_EQUAL_DEBUG( nCols, oldLocalData.getNumColumns() )
    DenseStorage<ValueType> newLocalData( rowDistribution->getLocalSize(), nCols );
    Redistributor redistributor( rowDistribution, getRowDistributionPtr() ); // target, source distributions
    redistributor.redistributeN( newLocalData.getData(), oldLocalData.getData(), nCols );
// COMMON_THROWEXCEPTION( "redistribution of dense rows not yet available" )
    oldLocalData.swap( newLocalData );
    this->setDistributionPtr( rowDistribution );
}

template<typename ValueType>
DenseMatrix<ValueType>::~DenseMatrix()
{
// Note: all member variables are freed by their own destructors
}

template<typename ValueType>
void DenseMatrix<ValueType>::setContextPtr( const ContextPtr context )
{
    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->setContextPtr( context );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getLocalRow( DenseVector<ValueType>& row, const IndexType iLocal ) const
{
    SCAI_LOG_INFO( logger, "get local row " << iLocal << " into " << row << " from this matrix: " << *this )
    SCAI_ASSERT_ERROR( row.getDistribution().isReplicated(), "row vector must be replicated" )
    const Distribution& distributionCol = getColDistribution();

    if ( distributionCol.isReplicated() )
    {
        // in this case we just can take the values from the local storage
        getLocalStorage().getRowImpl( row.getLocalValues(), iLocal );
        return;
    }

    // with column distribution: join the corresponding column data
    joinColumnData( row.getLocalValues(), iLocal, 1 );
    SCAI_LOG_INFO( logger, "joined ready" )
    SCAI_LOG_INFO( logger, "joined row value array = " << row.getLocalValues() )
    // just make sure that there is no mismatch of sizes
    SCAI_ASSERT_EQUAL_ERROR( row.getLocalValues().size(), row.size() );
}

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::getDiagonalImpl( DenseVector<OtherValueType>& diagonal ) const
{
    diagonal.allocate( getRowDistributionPtr() );
// const cast for local storage here is safe, otherwise we have to swap
    HArray<OtherValueType>& localValues = diagonal.getLocalValues();
    getLocalStorage().getDiagonal( localValues );
}

template<typename ValueType>
void DenseMatrix<ValueType>::getDiagonal( Vector& diagonal ) const
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

// todo: if ( diagonal.getVectorKind() == Vector::DENSE        )

    if ( true )
    {
// Dense vector with this row distribution, so we do not need a temporary array
        mepr::DenseMatrixWrapper<ValueType, SCAI_ARITHMETIC_HOST_LIST>::getDiagonalImpl( *this, diagonal );
        return;
    }

// Fallback solution with temporary arrays
    HArray<ValueType> localDiagonal;
    getLocalStorage().getDiagonal( localDiagonal );
    diagonal.assign( localDiagonal, getRowDistributionPtr() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const Vector& diagonal )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( getRowDistribution() != diagonal.getDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    getLocalStorage().setDiagonalV( diagonal.getLocalValues() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const Scalar diagonalValue )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    getLocalStorage().setDiagonal( diagonalValue.getValue<ValueType>() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::scale( const Vector& vector )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    getLocalStorage().scaleRows( vector.getLocalValues() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::scale( const Scalar scaleValue )
{
    getLocalStorage().scale( scaleValue.getValue<ValueType>() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::conj()
{
    getLocalStorage().conj();
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
Scalar DenseMatrix<ValueType>::getValue( IndexType i, IndexType j ) const
{
    ValueType myValue = static_cast<ValueType>( 0.0 );
    const Distribution& colDist = getColDistribution();
    const Distribution& rowDist = getRowDistribution();
    const Communicator& comm = rowDist.getCommunicator();

    if ( getRowDistribution().isLocal( i ) )
    {
        const IndexType iLocal = getRowDistribution().global2local( i );
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

        SCAI_ASSERT_ERROR( jLocal != nIndex, "non local column index" )
        SCAI_LOG_TRACE( logger,
                        "getting value for index(" << i << "," << j << ")" << " which is localy ( " << iLocal << "," << jLocal << " )" )
        myValue = mData[owner]->getValue( iLocal, jLocal );
    }

    SCAI_LOG_TRACE( logger, "My value is " << myValue << " starting sum reduction to produce final result." )
    return Scalar( comm.sum( myValue ) );
}

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesScalar( const Matrix& other, Scalar alpha )
{
    SCAI_LOG_INFO( logger, " this = " << alpha << " * " << other )
    assign( other );
    SCAI_LOG_INFO( logger, " this = other = " << *this )

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
    SCAI_REGION( "Mat.Dense.timesVector" )
    const HArray<ValueType>& localY = denseY.getLocalValues();
    HArray<ValueType>& localResult = denseResult.getLocalValues();
    ContextPtr localContext = mData[0]->getContextPtr();
    const Distribution& colDist = getColDistribution();
    const Communicator& comm = colDist.getCommunicator();
    int rank = comm.getRank();
    int n = colDist.getNumPartitions();
    mData[0]->prefetch();

    // It makes no sense to prefetch denseX because, if a transfer is started
    // the halo update needs to wait for this transfer to finish

    if ( betaValue != common::constants::ZERO )
    {
        denseY.prefetch( localContext );
    }

    const HArray<ValueType>& localX = denseX.getLocalValues();
    SCAI_LOG_INFO( logger,
                   comm << ": matrixTimesVector" << ", alpha = " << alphaValue << ", localX = " << localX << ", beta = " << betaValue << ", localY = " << localY )
    SCAI_LOG_INFO( logger,
                   "Aliasing: result = y : " << ( &denseResult == &denseY ) << ", local = " << ( &localResult == &localY ) )

    if ( n == 1 )
    {
// replicated column distribution, only on local block, X is replicated
// localResult = alpha * mData[0] * X + beta * localY
        const DenseStorage<ValueType>& dense = *mData[0];
        SCAI_LOG_INFO( logger, comm << ": matrixTimesVector, single dense block = " << dense )
        dense.matrixTimesVector( localResult, alphaValue, localX, betaValue, localY );
        return;
    }

    SCAI_LOG_INFO( logger, comm << ": start pipelined multiplication." )
    IndexType size = comm.max( localX.size() ); // largest local part of X
    mSendValues.clear();
    mReceiveValues.clear();
    ContextPtr contextPtr = Context::getHostPtr();
    HArray<ValueType>* sendValues = &mSendValues;
    HArray<ValueType>* recvValues = &mReceiveValues;
    {
// resize the receive buffer to be big enough for largest part of X
        WriteOnlyAccess<ValueType> wRecvValues( *recvValues, contextPtr, size );
        WriteOnlyAccess<ValueType> wSendValues( *sendValues, contextPtr, size );
        ReadAccess<ValueType> rLocalX( localX, contextPtr );
// fill send buffer with local X of this processor
        int i = 0;

        for ( ; i < localX.size(); ++i )
        {
            wSendValues[i] = rLocalX[i];
        }

        for ( ; i < size; ++i )
        {
            wSendValues[i] = static_cast<ValueType>( 0.0 );
        }
    }
    const int COMM_DIRECTION = 1; // shift buffer to next processor

    if ( Matrix::ASYNCHRONOUS == Matrix::getCommunicationKind() )
    {
        SCAI_LOG_INFO( logger, comm << ": asynchronous communication" )
// asynchronous communication always requires same sizes of arrays, might shift some more data
        common::unique_ptr<tasking::SyncToken> st( comm.shiftAsync( *recvValues, *sendValues, COMM_DIRECTION ) );
        SCAI_LOG_INFO( logger,
                       comm << ": matrixTimesVector, my dense block = " << *mData[rank] << ", localX = " << localX << ", localY = " << localY << ", localResult = " << localResult )
// overlap communication with local computation
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
                st.reset( comm.shiftAsync( *recvValues, *sendValues, COMM_DIRECTION ) );
            }
            else
            {
                st.reset( new tasking::NoSyncToken() );
            }

            SCAI_LOG_INFO( logger,
                           comm << ": matrixTimesVector, actual dense block [" << actualPartition << "] = " << *mData[actualPartition] << ", sendX = " << localX << ", localResult = " << localResult )
            // adapt the size of recvValues, that is now sendValues after swap
            HArray<ValueType> x( mData[actualPartition]->getNumColumns() );
            {
                static LAMAKernel<blaskernel::BLASKernelTrait::copy<ValueType> > copy;
                ContextPtr loc = this->getContextPtr();
                copy.getSupportedContext( loc );
                SCAI_CONTEXT_ACCESS( loc )
                ReadAccess<ValueType> readSend( *sendValues, loc );
                WriteAccess<ValueType> writeX( x, loc );
                copy[loc]( mData[actualPartition]->getNumColumns(), readSend.get(), 1, writeX.get(), 1 );
            }
            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, x, static_cast<ValueType>( 1.0 ), localResult );
            st->wait();
            std::swap( sendValues, recvValues );
        }
    }
    else
    {
// for synchronous communication we can use the real needed sizes
        {
            WriteAccess<ValueType> wSendValues( *sendValues );
            wSendValues.resize( localX.size() );
        }
        SCAI_LOG_INFO( logger, comm << ": synchronous communication" )
        comm.shiftArray( *recvValues, *sendValues, COMM_DIRECTION );
// For the synchronous shift we have no problems regarding the correct sizes
        SCAI_LOG_DEBUG( logger, comm << ": send " << *sendValues << ", recv " << *recvValues )
        SCAI_LOG_INFO( logger,
                       comm << ": matrixTimesVector, actual dense block [" << rank << "] = " << *mData[rank] << ", local X = " << localX << ", local Y = " << localY )
        mData[rank]->matrixTimesVector( localResult, alphaValue, localX, betaValue, localY );
        std::swap( sendValues, recvValues );

        for ( PartitionId p = 1; p < n; ++p )
        {
            PartitionId actualPartition = comm.getNeighbor( -p );
            comm.shiftArray( *recvValues, *sendValues, COMM_DIRECTION );
            SCAI_LOG_DEBUG( logger,
                            comm << ": send " << *sendValues << ", recv " << *recvValues << ", actual = " << actualPartition )
            SCAI_LOG_INFO( logger,
                           comm << ": matrixTimesVector, actual dense block [" << actualPartition << "] = " << *mData[actualPartition] << ", sendX = " << *sendValues << ", localResult = " << localResult )
            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, *sendValues, static_cast<ValueType>( 1.0 ), localResult );
            std::swap( sendValues, recvValues );
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::vectorTimesMatrixImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alphaValue,
    const DenseVector<ValueType>& denseX,
    const ValueType betaValue,
    const DenseVector<ValueType>& denseY ) const
{
    SCAI_REGION( "Mat.Dense.vectorTimesMatrix" )
    const HArray<ValueType>& localY = denseY.getLocalValues();
    HArray<ValueType>& localResult = denseResult.getLocalValues();
    ContextPtr localContext = mData[0]->getContextPtr();
    const Distribution& colDist = getColDistribution();
    const Communicator& comm = colDist.getCommunicator();
    mData[0]->prefetch();

    //It makes no sense to prefetch denseX because, if a transfer is started
    //the halo update needs to wait for this transfer to finish

    if ( betaValue != common::constants::ZERO )
    {
        denseY.prefetch( localContext );
    }

    const HArray<ValueType>& localX = denseX.getLocalValues();
    SCAI_LOG_INFO( logger,
                   comm << ": vectorTimesMatrix" << ", alpha = " << alphaValue << ", localX = " << localX << ", beta = " << betaValue << ", localY = " << localY )
    SCAI_LOG_INFO( logger,
                   "Aliasing: result = y : " << ( &denseResult == &denseY ) << ", local = " << ( &localResult == &localY ) )
    const DenseStorage<ValueType>& dense = *mData[0];
    SCAI_LOG_INFO( logger, comm << ": vectorTimesMatrix, singe dense block = " << dense )
    dense.vectorTimesMatrix( localResult, alphaValue, localX, betaValue, localY );
    return;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixPlusMatrix(
    const Scalar alpha,
    const Matrix& matA,
    const Scalar beta,
    const Matrix& matB )
{
    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", A = " << matA << ", B = " << matB )
    const DenseMatrix<ValueType>* denseA = dynamic_cast<const DenseMatrix<ValueType>*>( &matA );
    SCAI_ASSERT_ERROR( denseA, "Must be dense matrix<" << getValueType() << "> : " << matA )
    const DenseMatrix<ValueType>* denseB = dynamic_cast<const DenseMatrix<ValueType>*>( &matB );
    SCAI_ASSERT_ERROR( denseB, "Must be dense matrix<" << getValueType() << "> : " << matB )
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
    SCAI_REGION( "Mat.plusMatrix" )
// already verified
    SCAI_ASSERT_EQUAL_DEBUG( A.getRowDistribution(), B.getRowDistribution() )
    SCAI_ASSERT_EQUAL_DEBUG( A.getColDistribution(), B.getColDistribution() )
// Now we can do it completely local
    Matrix::setDistributedMatrix( A.getRowDistributionPtr(), A.getColDistributionPtr() );
// Add matrices of each chunk
    SCAI_LOG_INFO( logger, "Mat.plusMatrix, mDataSize = " << mData.size() );

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
    SCAI_ASSERT_ERROR( getRowDistribution().isReplicated(), "this->rows are distributed" )
    SCAI_ASSERT_ERROR( getColDistribution().isReplicated(), "this->cols are distributed" )
    SCAI_ASSERT_ERROR( B.getRowDistribution().isReplicated(), "B.rows are distributed" )
    SCAI_ASSERT_ERROR( B.getColDistribution().isReplicated(), "B.cols are distributed" )
    SCAI_ASSERT_ERROR( C.getRowDistribution().isReplicated(), "C.rows are distributed" )
    SCAI_ASSERT_ERROR( C.getColDistribution().isReplicated(), "C.cols are distributed" )
// Prefetch values to the ComputeLocation
    DenseMatrix* res = dynamic_cast<DenseMatrix*>( &result );

    if ( res == NULL )
    {
        COMMON_THROWEXCEPTION( "Only DenseMatrix DenseMatrix Multiplication is supported." )
    }

    const DenseMatrix* Bp = dynamic_cast<const DenseMatrix*>( &B );

    if ( Bp == NULL )
    {
        COMMON_THROWEXCEPTION( "Only DenseMatrix DenseMatrix Multiplication is supported." )
    }

    const DenseMatrix* Cp = dynamic_cast<const DenseMatrix*>( &C );

    if ( Cp == NULL )
    {
        COMMON_THROWEXCEPTION( "Only DenseMatrix DenseMatrix Multiplication is supported." )
    }

    if ( res == this )
    {
        SCAI_LOG_DEBUG( logger, "result is aliased with this A matrix" )
    }
    else if ( res == Bp )
    {
        SCAI_LOG_DEBUG( logger, "result is aliased with B matrix" )
    }
    else if ( res == Cp && beta.getValue<ValueType>() != 0.0 )
    {
        SCAI_LOG_DEBUG( logger, "result is aliased with C matrix" )
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "result is not aliased, so allocate it correctly" )
        res->allocate( getRowDistributionPtr(), B.getColDistributionPtr() );
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
    ValueType myMaxDiff = static_cast<ValueType>( 0.0 );

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        ValueType maxDiff = mData[i]->maxNorm();

        if ( maxDiff > myMaxDiff )
        {
            myMaxDiff = maxDiff;
        }
    }

    const Communicator& comm = getRowDistribution().getCommunicator();
    return Scalar( comm.max( myMaxDiff ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseMatrix<ValueType>::l1Norm() const
{
    const Communicator& comm = getRowDistribution().getCommunicator();
    ValueType mySum = static_cast<ValueType>( 0.0 );
    IndexType n = mData.size();

    for ( IndexType i = 0; i < n; i++ )
    {
        mySum += mData[i]->l1Norm();
    }

    return Scalar( comm.sum( mySum ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseMatrix<ValueType>::l2Norm() const
{
    const Communicator& comm = getRowDistribution().getCommunicator();
    ValueType mySum = static_cast<ValueType>( 0.0 );
    ValueType tmp;
    IndexType n = mData.size();

    for ( IndexType i = 0; i < n; i++ )
    {
        tmp = mData[i]->l2Norm();
        mySum += tmp * tmp;
    }

    return Scalar( common::Math::sqrt( comm.sum( mySum ) ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
Scalar DenseMatrix<ValueType>::maxDiffNorm( const Matrix& other ) const
{
    if ( !( ( mNumColumns == other.getNumColumns() ) && ( mNumRows == other.getNumRows() ) ) )
    {
        COMMON_THROWEXCEPTION( "maxDiffNorm requires matrices of same format" );
    }

// Implementation works only for same distributions and same type

    if ( ( getRowDistribution() == other.getRowDistribution() ) && ( getColDistribution() == other.getColDistribution() )
            && ( getValueType() == other.getValueType() ) )
    {
        const DenseMatrix<ValueType>* typedOther = dynamic_cast<const DenseMatrix<ValueType>*>( &other );
        SCAI_ASSERT_DEBUG( typedOther, "SERIOUS: wrong dynamic cast: " << other )
        return Scalar( maxDiffNormImpl( *typedOther ) );
    }
    else
    {
        SCAI_UNSUPPORTED( "maxDiffNorm requires temporary of " << other )
        DenseMatrix<ValueType> typedOther( other, getRowDistributionPtr(), getColDistributionPtr() );
        return Scalar( maxDiffNormImpl( typedOther ) );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseMatrix<ValueType>::maxDiffNormImpl( const DenseMatrix<ValueType>& other ) const
{
// implementation only supported for same distributions
    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution(), other.getRowDistribution() )
    SCAI_ASSERT_EQUAL_ERROR( getColDistribution(), other.getColDistribution() )
    ValueType myMaxDiff = static_cast<ValueType>( 0.0 );

    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        ValueType maxDiff = mData[i]->maxDiffNorm( *other.mData[i] );

        if ( maxDiff > myMaxDiff )
        {
            myMaxDiff = maxDiff;
        }
    }

    const Communicator& comm = getRowDistribution().getCommunicator();
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
void DenseMatrix<ValueType>::prefetch( hmemo::ContextPtr loc ) const
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
    SCAI_ASSERT_ERROR( mData.size() > 0, *this << ": no local values allocated" )

    if ( mData.size() == 1 )
    {
        return *mData[0];
    }

    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution(), getColDistribution() )
    const PartitionId myRank = getRowDistribution().getCommunicator().getRank();
    return *mData[myRank];
}

template<typename ValueType>
DenseStorage<ValueType>& DenseMatrix<ValueType>::getLocalStorage()
{
    SCAI_ASSERT_ERROR( mData.size() > 0, "no local values allocated" )

    if ( mData.size() == 1 )
    {
        return *mData[0];
    }

    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution(), getColDistribution() )
    const PartitionId myRank = getRowDistribution().getCommunicator().getRank();
    return *mData[myRank];
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getLocalNumValues() const
{
    // only locally stored number of values
    return getRowDistribution().getLocalSize() * mNumColumns;
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getLocalNumRows() const
{
    // only locally stored number of values
    return getRowDistribution().getLocalSize();
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

    return getRowDistribution().getCommunicator().sum( myNumValues );
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
    common::scalar::ScalarType type = common::getScalarType<ValueType>();
    stream << "DenseMatrix<" << type << ">( size = " << mNumRows << " x " << mNumColumns << ", rowdist = "
           << getRowDistribution() << ", coldist = " << getColDistribution() << ")";
}

template<typename ValueType>
common::scalar::ScalarType DenseMatrix<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

template<typename ValueType>
size_t DenseMatrix<ValueType>::getValueTypeSize() const
{
    return sizeof( ValueType );
}

template<typename ValueType>
size_t DenseMatrix<ValueType>::getMemoryUsage() const
{
    size_t memoryUsage = 0;

    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        memoryUsage += mData[i]->getMemoryUsage();
    }

    return getRowDistribution().getCommunicator().sum( memoryUsage );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>* DenseMatrix<ValueType>::newMatrix() const
{
    SCAI_LOG_INFO( logger, "SparseMatrix<ValueType>::newMatrix" )
    // use auto pointer for new sparse matrix to get data freed in case of Exception
    common::unique_ptr<DenseMatrix<ValueType> > newDenseMatrix( new DenseMatrix<ValueType>() );
    // inherit the context for local and halo storage
    newDenseMatrix->setContextPtr( this->getContextPtr() );
    newDenseMatrix->setCommunicationKind( this->getCommunicationKind() );
    SCAI_LOG_INFO( logger,
                   *this << ": create -> " << *newDenseMatrix << " @ " << * ( newDenseMatrix->getContextPtr() ) << ", kind = " << newDenseMatrix->getCommunicationKind() );
    return newDenseMatrix.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>* DenseMatrix<ValueType>::copy() const
{
    return new DenseMatrix<ValueType>( *this );
}

template<typename ValueType>
Matrix* DenseMatrix<ValueType>::create()
{
    return new DenseMatrix<ValueType>();
}

template<typename ValueType>
MatrixCreateKeyType DenseMatrix<ValueType>::createValue()
{
    common::scalar::ScalarType skind = common::getScalarType<ValueType>();
    return MatrixCreateKeyType ( Format::DENSE, skind );
}

template<typename ValueType>
MatrixCreateKeyType DenseMatrix<ValueType>::getCreateValue() const
{
    return createValue();
}

/* ========================================================================= */

template<typename ValueType>
const char* DenseMatrix<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
std::string DenseMatrix<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "DenseMatrix<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* DenseMatrix<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseMatrix, SCAI_ARITHMETIC_HOST )

} /* end namespace lama */

} /* end namespace scai */

