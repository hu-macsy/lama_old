/**
 * @file DenseMatrix.cpp
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
 * @brief DenseMatrix.cpp
 * @author Michael Drost
 * @date 22.02.2011
 */

// hpp
#include <scai/lama/matrix/DenseMatrix.hpp>

// local library
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/dmemo/NoDistribution.hpp>
#include <scai/dmemo/Redistributor.hpp>

// internal scai libraries
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/ScalarType.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::hmemo;
using namespace scai::dmemo;

namespace scai
{

using common::TypeTraits;
using utilskernel::LAMAKernel;

namespace lama
{

/* ========================================================================= */

template<typename ValueType, typename TList>
struct DenseMatrixWrapper;

template<typename ValueType>
struct DenseMatrixWrapper<ValueType, common::mepr::NullType>
{
    static void assignDense( DenseMatrix<ValueType>&, const _Matrix& other )
    {
        COMMON_THROWEXCEPTION( "type dense matrix not supported --> " << other )
    }
};

template<typename ValueType, typename H, typename T>
struct DenseMatrixWrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static void assignDense( DenseMatrix<ValueType>& obj, const _Matrix& other )
    {
        if ( other.getValueType() == common::getScalarType<H>() )
        {
            obj.copyDenseMatrix( reinterpret_cast<const DenseMatrix<H>& >( other ) );
        }
        else
        {
            DenseMatrixWrapper<ValueType, T>::assignDense( obj, other );
        }

    }
};

/* ========================================================================= */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseMatrix<ValueType>::logger, "Matrix.DenseMatrix" )

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::computeOwners()
{
    // enable any addressing for the column distribution

    getColDistribution().enableAnyAddressing();
}

/* ========================================================================= */
/*       Methods of DenseMatrix<ValueType>                                   */
/* ========================================================================= */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix()
{
    computeOwners();
    allocateData();      // will initialize it with zero
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const IndexType numRows, const IndexType numColumns ) :

    Matrix<ValueType>( numRows, numColumns )

{
    computeOwners();
    allocateData();      // will initialize it with zero
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DistributionPtr rowDist, DistributionPtr colDist ) :

    Matrix<ValueType>( rowDist, colDist )

{
    computeOwners();
    allocateData();      // will initialize it with zero
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
    const _Matrix& other,
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

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SMM_SM<ValueType>& expression )
{
    // resolve expression in base class matrix
    Matrix<ValueType>::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SMM<ValueType>& expression )
{
    // resolve expression in base class matrix
    Matrix<ValueType>::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SM_SM<ValueType>& expression )
{
    // resolve expression in base class matrix, usually -> matrixPlusMatrix
    Matrix<ValueType>::operator=( expression );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const Expression_SM<ValueType>& expression )
{
    // resolve expression in base class matrix
    Matrix<ValueType>::operator=( expression );
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
DenseMatrix<ValueType>::DenseMatrix( const _MatrixStorage& globalData )
{
    DistributionPtr rowDist( new NoDistribution( globalData.getNumRows() ) );
    DistributionPtr colDist( new NoDistribution( globalData.getNumColumns() ) );
    DenseMatrix<ValueType>::assign( globalData, rowDist, colDist );
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
DenseMatrix<ValueType>::DenseMatrix( const DenseMatrix<ValueType>& other ) :

    Matrix<ValueType>()

{
    SCAI_LOG_INFO( logger, "copy constructor( dense matrix, same value type) : " << other )
    assign( other ); // will choose the local assignment
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const _Matrix& other, bool transposeFlag )
{
    SCAI_LOG_INFO( logger, "copy constructor( any matrix) : " << other << ", transpse = " << transposeFlag )

    if ( transposeFlag )
    {
        assignTranspose( other );
    }
    else
    {
        assign( other );
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setIdentity( DistributionPtr dist )
{
    _Matrix::setDistributedMatrix( dist, dist );
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

    _Matrix::setDistributedMatrix( rowDist, tmpReplicatedColDistribution );
    // due to temporary replicated col distribution, mData has only one entry
    mData[0]->setDenseData( n, m, values, eps.getValue<ValueType>() );
    SCAI_LOG_INFO( logger,
                   "Dense matrix, row dist = " << *rowDist << " filled locally with " << ( n * m ) << " values, now split for col dist = " << *colDist );

    if ( !colDist->isReplicated() )
    {
        splitColumns( colDist );

        for ( PartitionId i = 0; i < colDist->getCommunicator().getSize(); ++i )
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

    _Matrix::setDistributedMatrix( rowDist, tmpReplicatedColDistribution );
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
    DistributionPtr rowDist,
    DistributionPtr colDist,
    const IndexType numDiagonals,
    const HArray<IndexType>& offsets,
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

    _Matrix::setDistributedMatrix( rowDist, tmpReplicatedColDistribution );

    // due to temporary replicated col distribution, mData has only one entry

    mData[0]->setDIAData( n, m, numDiagonals, offsets, values );

    if ( !colDist->isReplicated() )
    {
        splitColumns( colDist );
    }
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
        _Matrix::checkSettings();

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
void DenseMatrix<ValueType>::invert( const _Matrix& other )
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
    mData[0]->setCSRData( rowIA.size() + 1, getNumColumns(), rowJA.size(), rowIA, rowJA, rowValues );
    // ToDo: split up mData[0] according to column distribution
}

template<typename ValueType>
void DenseMatrix<ValueType>::clear()
{
    _Matrix::setReplicatedMatrix( 0, 0 ); // clear _Matrix
    mData.resize( 1 ); // clear Data
    mData[0]->clear();
}

template<typename ValueType>
void DenseMatrix<ValueType>::purge()
{
    _Matrix::setReplicatedMatrix( 0, 0 ); // clear _Matrix
    mData.resize( 1 ); // clear Data
    mData[0]->purge();
}

template<typename ValueType>
void DenseMatrix<ValueType>::allocate( const IndexType numRows, const IndexType numColumns )
{
    DistributionPtr rowDist( new NoDistribution( numRows ) );
    DistributionPtr colDist( new NoDistribution( numColumns ) );

    allocate( rowDist, colDist );
}

template<typename ValueType>
void DenseMatrix<ValueType>::allocate( DistributionPtr rowDistribution, DistributionPtr colDistribution )
{
    SCAI_LOG_INFO( logger,
                   *this << " with mData[" << mData.size() << "]" << ", allocate row dist = " << *rowDistribution 
                         << ", col dist = " << *colDistribution )

    _Matrix::setDistributedMatrix( rowDistribution, colDistribution );

    computeOwners();
    allocateData();

    SCAI_LOG_DEBUG( logger, *this << ": now allocated" )
}

template<typename ValueType>
void DenseMatrix<ValueType>::swap( DenseMatrix<ValueType>& other )
{
    _Matrix::swapMatrix( other );
    // now swap own member variables
    std::swap( mData, other.mData );
}

template<typename ValueType>
void DenseMatrix<ValueType>::assignTranspose( const _Matrix& other  )
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
    SCAI_REGION( "Mat.Dense.assignTranspose" )

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
void DenseMatrix<ValueType>::assign( const _Matrix& other )
{
    SCAI_LOG_INFO( logger, "assign " << other << " to " << *this )

    if ( &other == this )
    {
        SCAI_LOG_INFO( logger, "self assign, is skpped" )
    }
    else if ( other.getMatrixKind() == MatrixKind::DENSE )
    {
        SCAI_LOG_INFO( logger, "copy dense matrix" )
        DenseMatrixWrapper<ValueType, SCAI_NUMERIC_TYPES_HOST_LIST>::assignDense( *this, other );
    }
    else if ( other.getMatrixKind() == MatrixKind::SPARSE )
    {
        SCAI_LOG_INFO( logger, "copy sparse matrix" )
        assignSparse( other );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Unsupported: assign " << other << " to " << *this )
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignSparse( const _Matrix& other )
{
    // we need replicated column distribution to get this routine working

    if ( !other.getColDistribution().isReplicated() )
    {
        DistributionPtr repColDist( new NoDistribution( other.getNumColumns() ) );

        // std::unique_ptr<Matrix> tmpOther( other.copy() );
        // tmpOther->redistribute( other.getRowDistributionPtr(), repColDist );
        // SCAI_LOG_WARN( logger, "create temporary matrix with replicated columns: " << *tmpOther )
        // assignSparse( *tmpOther );

        CSRSparseMatrix<ValueType> repOther( other, other.getRowDistributionPtr(), repColDist );

        assignSparse( repOther );

        splitColumns( other.getColDistributionPtr() );

        return;
    }

    // replicated columns in sparse matrix, so we can assign local data

    _Matrix::setDistributedMatrix( other.getRowDistributionPtr(), other.getColDistributionPtr() );
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
    _Matrix::setReplicatedMatrix( numRows, numColumns );
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( storage ) );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist )
{
    SCAI_LOG_INFO( logger, "assign local matrix storage = " << storage )

    _Matrix::setDistributedMatrix( rowDist, colDist );

    colDist->enableAnyAddressing();

    if ( storage.getNumRows() == rowDist->getLocalSize() )
    {
        // only format conversion of the local storage, @todo avoid it if storage is DenseStorage<ValueType>

        if ( storage.getFormat() == Format::DENSE && storage.getValueType() == getValueType() )
        {
            const DenseStorage<ValueType>* localData = dynamic_cast<const DenseStorage<ValueType>*>( &storage );
            SCAI_ASSERT_ERROR( localData, "dynamic_cast<constDenseStorage<ValueType>*> failed: " << storage )

            splitColumnData( mData, *localData, *colDist );
        }
        else if ( colDist->isReplicated() )
        {
            std::shared_ptr<DenseStorage<ValueType> > dataPtr( new DenseStorage<ValueType>( storage ) );

            mData.resize( 1 );
            mData[0] = dataPtr;
        }
        else
        {
            DenseStorage<ValueType> localData( storage );
            splitColumnData( mData, localData, *colDist );
        }
    }
    else if ( storage.getNumRows() == rowDist->getGlobalSize() )
    {
        // we also localize the rows of the matrix

        DenseStorage<ValueType> localData;
        localData.localize( storage, *rowDist );
        splitColumnData( mData, localData, *colDist );
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
    SCAI_LOG_INFO( logger, "build local storage with replicated columns for " << *this )

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
        DenseStorage<ValueType> denseStorage( numLocalRows, getNumColumns() );
        joinColumnData( denseStorage.getData(), 0, numLocalRows );
        storage = denseStorage;
    }

    SCAI_LOG_DEBUG( logger, "buildLocalStorage( " << *this << " ) = " << storage )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::joinColumnData(
    HArray<ValueType>& result,
    const IndexType firstRow,
    const IndexType nRows ) const
{
    SCAI_REGION( "Mat.Dense.joinColumnData" )

    SCAI_LOG_INFO( logger, "join column data, firstRow = " << firstRow << ", nRows = " << nRows )

    const IndexType nColumns = getNumColumns();

    const PartitionId numColPartitions = static_cast<PartitionId>( mData.size() );

    ContextPtr hostContext = Context::getHostPtr();

    HArray<IndexType> offsets;       // running sizes for each partition
    HArray<IndexType> local2global;  // global indexes sorted by owners

    getColDistribution().getAnyLocal2Global( offsets, local2global );

    SCAI_ASSERT_EQ_DEBUG( offsets.size(), numColPartitions + 1, "serious mismatch" )
    SCAI_ASSERT_EQ_DEBUG( local2global.size(), nColumns, "serious mismatch" )

    {
        WriteOnlyAccess<ValueType> wResult( result, hostContext, nRows * nColumns );

        ReadAccess<IndexType> rLocal2Global( local2global, hostContext );
        ReadAccess<IndexType> rOffsets( offsets, hostContext );

        // scatter local data of each chunk into the full result array

        for ( PartitionId p = 0; p < numColPartitions; ++p )
        {
            ReadAccess<ValueType> rLocalData( mData[p]->getData(), hostContext );

            IndexType offset = rOffsets[p];
            IndexType numLocalColumns = rOffsets[p + 1] - rOffsets[p];

            for ( IndexType localJ = 0; localJ < numLocalColumns; ++localJ )
            {
                IndexType globalJ = rLocal2Global[ offset + localJ ];
    
                for ( IndexType i = firstRow; i < firstRow + nRows ; ++i )
                {
                    wResult[ globalJ + ( i - firstRow ) * nColumns] = rLocalData[ localJ + i * numLocalColumns ];
                }
            }
        }
    }
 
    SCAI_LOG_DEBUG( logger, "ready join column data" )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::allocateData()
{
    const Distribution& colDist = getColDistribution();

    SCAI_ASSERT_EQUAL_DEBUG( getNumColumns(), colDist.getGlobalSize() )

    const PartitionId numChunks = colDist.getCommunicator().getSize();

    mData.resize( numChunks );

    const IndexType numLocalRows = getRowDistribution().getLocalSize();

    SCAI_LOG_INFO( logger, "build " << numChunks << " data arrays for numLocalRows = " << numLocalRows );

    if ( numChunks == 1 )
    {
        if ( mData[0] )
        {
            // just reallocate the storage
            mData[0]->allocate( numLocalRows, getNumColumns() );
        }
        else
        {
            // first time allocation
            mData[0].reset( new DenseStorage<ValueType>( numLocalRows, getNumColumns() ) );
        }

        return;
    }

    ContextPtr ctx = Context::getHostPtr();

    IndexType count = 0;   // sum up the sizes, verify correct sum

    for ( PartitionId p = 0; p < numChunks; ++p )
    {
        IndexType numLocalColumns = colDist.getAnyLocalSize( p );
        count += numLocalColumns;
        mData[p].reset( new DenseStorage<ValueType>( numLocalRows, numLocalColumns ) );
    }

    SCAI_ASSERT_EQ_ERROR( count, getNumColumns(), "Illegal owners." )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::splitColumnData(
    std::vector<std::shared_ptr<DenseStorage<ValueType> > >& chunks,
    const DenseStorage<ValueType>& columnData,
    const Distribution& colDistribution )
{
    SCAI_ASSERT_EQ_ERROR( colDistribution.getGlobalSize(), columnData.getNumColumns(), "serious mismatch" )

    if ( colDistribution.isReplicated() )
    {
        SCAI_REGION( "Mat.Dense.noSplitColumnData" )

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

    SCAI_REGION( "Mat.Dense.splitColumnData" )

    PartitionId numPartitions = colDistribution.getNumPartitions();

    SCAI_LOG_INFO( logger, "split columns of " << columnData << " into " << numPartitions << " chunks" )

    const IndexType numColumns = columnData.getNumColumns();
    const IndexType numRows = columnData.getNumRows();

    HArray<IndexType> offsets;       // running sizes for each partition
    HArray<IndexType> local2global;  // global indexes sorted by owners

    colDistribution.getAnyLocal2Global( offsets, local2global );

    SCAI_ASSERT_EQ_DEBUG( local2global.size(), numColumns, "serious mismatch" )
    SCAI_ASSERT_EQ_DEBUG( offsets.size(), numPartitions + 1, "serious mismatch" )

    chunks.clear();
    chunks.resize( numPartitions );

    ContextPtr ctx = Context::getHostPtr();

    ReadAccess<ValueType> columnDataRead( columnData.getData(), ctx );
    ReadAccess<IndexType> rLocal2Global( local2global, ctx );
    ReadAccess<IndexType> rOffsets( offsets, ctx );

    // gather the data for each partition from the global data

    for ( PartitionId p = 0; p < numPartitions; ++p )
    {
        IndexType offset = rOffsets[p];
        IndexType numLocalColumns = rOffsets[p + 1] - offset;

        chunks[p].reset( new DenseStorage<ValueType>( numRows, numLocalColumns ) );

        WriteAccess<ValueType> wChunkData( chunks[p]->getData(), ctx );

        IndexType pos = 0;  // traversing the elements of chunk data for p-th partition

        for ( IndexType i = 0; i < numRows; ++i )
        {
            for ( IndexType localCol = 0; localCol < numLocalColumns; ++localCol )
            {
                IndexType globalCol = rLocal2Global[localCol + offset];
                wChunkData[pos++] = columnDataRead[ i * numColumns + globalCol ];
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

    // As only rows are exchanged, col distribution must be replicated

    if ( getColDistribution().getNumPartitions() != 1 )
    {
        // Join all column data
        const IndexType numCols = getNumColumns();
        const IndexType numLocalRows = getRowDistribution().getLocalSize();
        std::shared_ptr<DenseStorage<ValueType> > colData;
        colData.reset( new DenseStorage<ValueType>( numLocalRows, numCols ) );
        joinColumnData( colData->getData(), 0, numLocalRows );
        mData.clear();
        mData.resize( 1 );
        mData[0] = colData;
        DistributionPtr noColDist( new NoDistribution( getNumColumns() ) );
        _Matrix::setDistributedMatrix( getRowDistributionPtr(), noColDist );
    }

    redistributeRows( rowDistribution );
    splitColumns( colDistribution );
}

template<typename ValueType>
void DenseMatrix<ValueType>::splitColumns( DistributionPtr colDistribution )
{
    SCAI_ASSERT_EQUAL_ERROR( 1, getColDistribution().getNumPartitions() )
    std::shared_ptr<DenseStorage<ValueType> > oldStorage = mData[0];
    _Matrix::setDistributedMatrix( getRowDistributionPtr(), colDistribution );
    computeOwners(); // compute mapping column index -> chunk
    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution().getLocalSize(), oldStorage->getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), oldStorage->getNumColumns() )
    splitColumnData( mData, *oldStorage, *colDistribution );
// old storage will be freed here at end of scope
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::redistribute( const Redistributor& redistributor, DistributionPtr colDistributionPtr )
{
    SCAI_ASSERT_EQ_ERROR( getRowDistribution(), *redistributor.getSourceDistributionPtr(),
                          "redistributor does not match to actual distribution of this vector" );

    if ( !getColDistribution().isReplicated() )
    {
        // Halo must be removed before redistribution 

        DistributionPtr repColDistributionPtr( new NoDistribution( getNumColumns() ) );
        redistribute( getRowDistributionPtr(), repColDistributionPtr );
    }

    std::shared_ptr<DenseStorage<ValueType> > newData( mData[0]->newMatrixStorage() );
    newData->redistribute( *mData[0], redistributor );
    mData[0] = newData;

    _Matrix::setDistributedMatrix( redistributor.getTargetDistributionPtr(), getColDistributionPtr() );

    redistribute( getRowDistributionPtr(), colDistributionPtr );
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
        DenseStorage<ValueType> newLocalData( getNumRows(), nCols );
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

/* -------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::~DenseMatrix()
{
    // nothing to do, all member variables are freed by their own destructors
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::setContextPtr( const ContextPtr context )
{
    SCAI_ASSERT_ERROR( context.get(), "NULL context for dense matrix, not allowed" )

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->setContextPtr( context );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getLocalRow( HArray<ValueType>& row, const IndexType localRowIndex ) const
{
    SCAI_REGION( "Mat.Dense.getLocalRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( localRowIndex, getRowDistribution().getLocalSize(), "illegal local row index" );

    SCAI_LOG_INFO( logger, "get local row " << localRowIndex << " from this matrix: " << *this )

    const Distribution& distributionCol = getColDistribution();

    if ( distributionCol.isReplicated() )
    {
        // in this case we just can take the values from the local storage
        getLocalStorage().getRowImpl( row, localRowIndex );
        return;
    }

    // with column distribution: join the corresponding column data

    joinColumnData( row, localRowIndex, 1 );

    SCAI_LOG_INFO( logger, "local row with joined column data: " << row )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getRowLocal( _Vector&, const IndexType ) const
{
    COMMON_THROWEXCEPTION( "not available yet" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getRow( _Vector& row, const IndexType globalRowIndex ) const
{
    // if v is not a dense vector or not of same type, use a temporary dense vector

    if ( row.getVectorKind() != VectorKind::DENSE || row.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "getRow requires temporary" )
        DenseVector<ValueType> denseRow;
        getRow( denseRow, globalRowIndex );
        row.assign( denseRow );   // transform the dense vector into sparse vector
        return;
    }

    SCAI_REGION( "Mat.Dense.getRow" )

    DenseVector<ValueType>& denseRow = reinterpret_cast<DenseVector<ValueType>&>( row );

    denseRow.allocate( getColDistributionPtr() );   // same dist as column dist

    HArray<ValueType>& values = denseRow.getLocalValues();  // be careful to guarantee consistency

    if ( getRowDistribution().isReplicated() )
    {
        SCAI_LOG_INFO( logger, "getRow with replicated row distribution" )

        // each processor has all data, just pick up my local part

        const Communicator& comm = getColDistribution().getCommunicator();

        PartitionId colRank = comm.getRank();

        mData[ colRank ]->getRow( values, globalRowIndex );

        SCAI_ASSERT_EQ_ERROR( values.size(), getColDistribution().getLocalSize(), "serious mismatch" )

        return;
    }

    const Communicator& comm = getRowDistribution().getCommunicator();

    // Note: for the row distribution any owner might not be enabled

    PartitionId rowOwner = getRowDistribution().findOwner( globalRowIndex );

    SCAI_LOG_INFO( logger, "row dist = " << getRowDistribution() 
                            << ", owner = " << rowOwner << " for row " << globalRowIndex )

    PartitionId np = getColDistribution().getNumPartitions();

    if ( np == 1 )
    {
        SCAI_LOG_DEBUG( logger, "getRow with replicated col distribution, row owner = " << rowOwner )

        // owner gets the dense row locally and broadcasts it

        if ( rowOwner == comm.getRank() )
        {
            IndexType localRowIndex = getRowDistribution().global2local( globalRowIndex );
            mData[0]->getRow( values, localRowIndex );
        }

        getRowDistribution().getCommunicator().bcastArray( values, rowOwner );

        return;
    }

    SCAI_LOG_DEBUG( logger, comm << ": bcast chunks for this col dist = " << getColDistribution() )

    if ( rowOwner == comm.getRank() )
    {
        IndexType localRowIndex = getRowDistribution().global2local( globalRowIndex );

        SCAI_ASSERT_EQ_ERROR( static_cast<IndexType>( mData.size() ), np, "illegal column data" )

        HArray<ValueType> sendBuffer;
        HArray<ValueType> recvBuffer;

        CommunicationPlan sendPlan;
        CommunicationPlan recvPlan( NULL, 0 );  // nothing to receive

        for ( PartitionId p = 0; p < comm.getSize(); ++p )
        {
            if ( p == comm.getRank() )
            {
                // pick up the local values

                mData[p]->getRow( values, localRowIndex );
            }
            else
            {
                mData[p]->getRow( sendBuffer, localRowIndex );
             
                if ( sendBuffer.size() )
                {
                    sendPlan.singleEntry( p, sendBuffer.size() );
                    comm.exchangeByPlan( recvBuffer, recvPlan, sendBuffer, sendPlan );
                }
            }
        }
    }
    else
    {
        HArray<ValueType> dummySend;

        // Not owner, build recv plan

        IndexType size = getColDistribution().getLocalSize();

        CommunicationPlan sendPlan( NULL, 0 );
        CommunicationPlan recvPlan( NULL, 0 );

        recvPlan.singleEntry( rowOwner, size );

        SCAI_LOG_DEBUG( logger, comm << ": getRow, recvPlan = " << recvPlan << ", sendPlan = " << sendPlan )

        comm.exchangeByPlan( values, recvPlan, dummySend, sendPlan );
    }

    // guarantee consistency in the dense vector for the local data

    SCAI_ASSERT_EQ_ERROR( values.size(), getColDistribution().getLocalSize(), "serious mismatch" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getColumn( _Vector& col, const IndexType globalColIndex ) const
{
    // if col is not a dense vector, use a temporary dense vector

    if ( col.getVectorKind() != VectorKind::DENSE || col.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "getCol requires temporary, use DenseVector on DenseMatrix" )
        DenseVector<ValueType> denseColumn;
        getColumn( denseColumn, globalColIndex );
        col.assign( denseColumn );   // transform the dense vector into sparse vector, works for all
        return;
    }

    SCAI_REGION( "Mat.Dense.getColumn" )

    SCAI_ASSERT_DEBUG( dynamic_cast<DenseVector<ValueType>*>( &col ), "col not DenseVector<" << getValueType() << ">" )

    DenseVector<ValueType>& denseCol = reinterpret_cast<DenseVector<ValueType>&>( col );

    // result vector inherits the row distribution 

    denseCol.allocate( getRowDistributionPtr() );

    // find the owner and local column index of col

    PartitionId owner         = 0;
    IndexType   localColIndex = globalColIndex;

    const Distribution& colDist = getColDistribution();

    if ( !colDist.isReplicated() )
    {
        // determine the owner and local index 

        owner = colDist.getAnyOwner( globalColIndex );
        localColIndex = colDist.getAnyLocalIndex( globalColIndex, owner );
    }

    SCAI_ASSERT_DEBUG( mData[owner], "No data for owner = " << owner )

    SCAI_LOG_INFO( logger, "getColumn( " << globalColIndex << " ) : owner = " << owner
                           << ", local col = " << localColIndex << ", mData = " << *mData[owner] )

    HArray<ValueType>& values = denseCol.getLocalValues();

    mData[owner]->getColumn( values, localColIndex );

    // verify that local size of data matches local size of distribution, so we have consistency

    SCAI_ASSERT_EQ_DEBUG( values.size(), denseCol.getDistribution().getLocalSize(), "serious mismatch" );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::setLocalRow(
    const hmemo::HArray<ValueType>& row,
    const IndexType localRowIndex,
    const common::BinaryOp op )
{
    SCAI_REGION( "Mat.Dense.setLocalRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( localRowIndex, getRowDistribution().getLocalSize(), "illegal local row index" )
    SCAI_ASSERT_EQ_DEBUG( row.size(), getNumColumns(), "size of row illegal" )

    const PartitionId numColPartitions = static_cast<PartitionId>( mData.size() );

    if ( numColPartitions == 1 )
    {
        // in this case we just can take the values from the local storage

        mData[0]->setRowImpl( row, localRowIndex, op );
        return;
    }

    utilskernel::LArray<IndexType> offsets;
    utilskernel::LArray<IndexType> perm;

    getColDistribution().getAnyLocal2Global( offsets, perm );

    HArray<ValueType> rowResorted;   // row resorted according to the owners

    utilskernel::HArrayUtils::gatherImpl( rowResorted, row, perm, common::BinaryOp::COPY );

    ReadAccess<IndexType> rOffsets( offsets );

    // now sort in the different parts

    for ( PartitionId ip = 0; ip < numColPartitions; ++ip )
    {
        // tricky workaround for: HArraySection<ValueType>( rowResorted, offset = .., inc = 1, n = ... )

        const ValueType* ptrRowPart;

        IndexType nRowPart = rOffsets[ ip + 1 ] - rOffsets[ip];

        {
            ReadAccess<ValueType> rRow( rowResorted );
            ptrRowPart = rRow.get() + rOffsets[ip];
            nRowPart   = rOffsets[ ip + 1 ] - rOffsets[ ip ];
        }

        HArrayRef<ValueType> rowPartition( nRowPart, ptrRowPart );

        mData[ip]->setRowImpl( rowPartition, localRowIndex, op );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::setLocalColumn(
    const hmemo::HArray<ValueType>& column,
    const IndexType globalColIndex,
    const common::BinaryOp op )
{
    SCAI_REGION( "Mat.Dense.setLocalColumn" )

    // find the owner and local column index of col

    PartitionId owner         = 0;
    IndexType   localColIndex = globalColIndex;

    const Distribution& colDist = getColDistribution();

    if ( !colDist.isReplicated() )
    {
        owner = colDist.getAnyOwner( globalColIndex );
        localColIndex = colDist.getAnyLocalIndex( globalColIndex, owner );
    }

    SCAI_ASSERT_ERROR( mData[owner], "No data for owner = " << owner )

    SCAI_LOG_INFO( logger, "setLocalColumn( " << globalColIndex << " ) : owner = " << owner
                   << ", local col = " << localColIndex  << ", mData = " << *mData[owner] )

    mData[owner]->setColumnImpl( column, localColIndex, op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getDiagonal( _Vector& diagonal ) const
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( diagonal.getVectorKind() != VectorKind::DENSE || diagonal.getValueType() != getValueType() )
    {
        DenseVector<ValueType> tmpDiagonal( diagonal.getContextPtr() );
        getDiagonal( tmpDiagonal );
        diagonal.assign( tmpDiagonal );   // does the correct type / kind conversion
        return;
    }

    // we can recast it now to dense vector, so we have access to its local values

    DenseVector<ValueType>& denseDiagonal = reinterpret_cast<DenseVector<ValueType>&>( diagonal );

    denseDiagonal.allocate( getRowDistributionPtr() );
    getLocalStorage().getDiagonal( denseDiagonal.getLocalValues() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const _Vector& diagonal )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "setDiagonal only for square matrices with same row/col distribution" )
    }

    if ( getRowDistribution() != diagonal.getDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( diagonal.getVectorKind() != VectorKind::DENSE || diagonal.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "setDiagonal: diagonal will be converted" )
        DenseVector<ValueType> tmpDiagonal( diagonal );
        setDiagonal( tmpDiagonal );
        return;
    }

    const DenseVector<ValueType>& diagonalDense = reinterpret_cast<const DenseVector<ValueType>&>( diagonal );

    getLocalStorage().setDiagonalV( diagonalDense.getLocalValues() );
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

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::reduce(
    _Vector& v, 
    const IndexType dim, 
    const common::BinaryOp reduceOp, 
    const common::UnaryOp elemOp ) const
{
    SCAI_REGION( "Mat.Dense.reduce" )

    // SCAI_ASSERT_EQ_ERROR( v.getValueType(), 

    DenseVector<ValueType>& denseV = reinterpret_cast<DenseVector<ValueType>&>( v );

    if ( dim == 0 )
    {
        denseV.allocate( getRowDistributionPtr() );

        denseV = ValueType( 0 );   // initialize v with neutral element

        for ( size_t k = 0; k < mData.size(); ++k )
        {
            mData[k]->reduce( denseV.getLocalValues(), 0, reduceOp, elemOp );
        }

        return;
    }

    if ( dim == 1 )
    {
        denseV.allocate( getColDistributionPtr() );

        denseV = ValueType( 0 );   // initialize v with neutral element

        if ( getRowDistribution().getCommunicator().getSize() == 1 )
        {
            // full matrix is replicated, the columns might have any distribution

            PartitionId rank = getColDistribution().getCommunicator().getRank();

            mData[rank]->reduce( denseV.getLocalValues(), 1, reduceOp, elemOp );

            return;   // matrix is replicated, compute just my values
        }

        // rows are distributed

        IndexType np = getColDistribution().getCommunicator().getSize();

        if ( np == 1 )
        {
             SCAI_ASSERT_EQ_ERROR( reduceOp, common::BinaryOp::ADD, "only add supported" )

             mData[0]->reduce( denseV.getLocalValues(), 1, reduceOp, elemOp );
             getRowDistribution().getCommunicator().sumArray( denseV.getLocalValues() );
             return;
        }

        // circular shift of local parts and reduce vals from here

        const int COMM_DIRECTION = 1;       // circular shifting
        HArray<ValueType>& sendValues = mSendValues;
        HArray<ValueType>& recvValues = mReceiveValues;
        IndexType maxSize = getColDistribution().getMaxLocalSize();

        // send/recv buffers must be large enough to keep largest amout of data from any processor

        ContextPtr contextPtr = Context::getHostPtr();
        recvValues.reserve( contextPtr, maxSize );
        sendValues.reserve( contextPtr, maxSize );

        utilskernel::HArrayUtils::assign( sendValues, denseV.getLocalValues() );

        const Communicator& comm = getColDistribution().getCommunicator();

        for ( PartitionId p = 0; p < np; ++p )
        {
            // compute the owner of the values that are in the current send buffer
            PartitionId actualPartition = comm.getNeighbor( -p );
            mData[actualPartition]->reduce( sendValues, 1, reduceOp, elemOp );
            comm.shiftArray( recvValues, sendValues, COMM_DIRECTION );
            std::swap( sendValues, recvValues );
        }

        utilskernel::HArrayUtils::assign( denseV.getLocalValues(), sendValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "illegal reduce dim = " << dim << " for dense matrix" )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::scale( const _Vector& vector )
{
    if ( getRowDistribution() != vector.getDistribution() )
    {
        COMMON_THROWEXCEPTION( "scale vector must have same distribution as matrix row distribution" )
    }

    if ( vector.getVectorKind() != VectorKind::DENSE || vector.getValueType() != getValueType() )
    {
        SCAI_LOG_WARN( logger, "scale: vector requires temporary" )
        DenseVector<ValueType> tmpVector( vector );
        scale( tmpVector );
        return;
    }
    
    const DenseVector<ValueType>& denseVector = reinterpret_cast<const DenseVector<ValueType>&>( vector );

    getLocalStorage().scaleRows( denseVector.getLocalValues() );
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
ValueType DenseMatrix<ValueType>::getValue( IndexType i, IndexType j ) const
{
    ValueType myValue = 0;
    const Distribution& colDist = getColDistribution();
    const Distribution& rowDist = getRowDistribution();
    const Communicator& commRow = rowDist.getCommunicator();

    if ( getRowDistribution().isLocal( i ) )
    {
        const IndexType iLocal = getRowDistribution().global2local( i );

        PartitionId owner = 0;
        IndexType  jLocal = nIndex;

        if ( colDist.getNumPartitions() == 1 )
        {
            owner  = 0;
            jLocal = j;
        }
        else
        {
            owner  = colDist.getAnyOwner( j );
            jLocal = colDist.getAnyLocalIndex( j, owner );
        }

        SCAI_ASSERT_ERROR( jLocal != nIndex, "non local column index" )
        SCAI_LOG_TRACE( logger,
                        "getting value for index(" << i << "," << j << ")" << " which is localy ( " << iLocal << "," << jLocal << " )" )
        myValue = mData[owner]->getValue( iLocal, jLocal );
    }

    SCAI_LOG_TRACE( logger, "My value is " << myValue << " starting sum reduction to produce final result." )

    return commRow.sum( myValue );
}

template<typename ValueType>
void DenseMatrix<ValueType>::setValue(
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

    PartitionId owner  = 0;

    IndexType   jLocal = nIndex;

    if ( distributionCol.getNumPartitions() == 1 )
    {
        jLocal = j;
    }
    else 
    {
        owner  = distributionCol.getAnyOwner( j );
        jLocal = distributionCol.getAnyLocalIndex( j, owner );
    }

    SCAI_ASSERT_ERROR( jLocal != nIndex, "non local column index" )

    mData[owner]->setValue( iLocal, jLocal, val, op );
}

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesScalar( const _Matrix& other, Scalar alpha )
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
    PartitionId rank = comm.getRank();
    PartitionId n = colDist.getNumPartitions();
    mData[0]->prefetch();

    // It makes no sense to prefetch denseX because, if a transfer is started
    // the halo update needs to wait for this transfer to finish

    if ( betaValue != common::Constants::ZERO )
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
        IndexType i = 0;

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

    if ( SyncKind::ASYNCHRONOUS == _Matrix::getCommunicationKind() )
    {
        SCAI_LOG_INFO( logger, comm << ": asynchronous communication" )
// asynchronous communication always requires same sizes of arrays, might shift some more data
        std::unique_ptr<tasking::SyncToken> st( comm.shiftAsync( *recvValues, *sendValues, COMM_DIRECTION ) );
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
            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, *sendValues, ValueType( 1 ), localResult );
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

    const Distribution& colDist = getColDistribution();
    const Distribution& rowDist = getRowDistribution();

    const Communicator& comm = colDist.getCommunicator();

    const HArray<ValueType>& localX = denseX.getLocalValues();

    PartitionId nParts = colDist.getNumPartitions();

    SCAI_ASSERT_GT( nParts, 1, "replicated columns are not considered here" )

    if ( rowDist.getNumPartitions() == 1 )
    {
        // the full matrix is replicated, the result is distributed, compute just its part

        PartitionId rank = comm.getRank();
        mData[rank]->vectorTimesMatrix( localResult, alphaValue, localX, betaValue, localY );

        SCAI_LOG_INFO( logger, comm << ": computed local Result for replicated matrix: " << localResult );

        return;
    }

    const int COMM_DIRECTION = 1;       // circular shifting

    // reuse member variables of this DenseMatrix, avoids too much reallocation

    HArray<ValueType>& sendValues = mSendValues;

    HArray<ValueType>& recvValues = mReceiveValues;

    IndexType maxSize = colDist.getMaxLocalSize();

    // send/recv buffers must be large enough to keep largest amout of data from any processor

    ContextPtr contextPtr = Context::getHostPtr();

    recvValues.reserve( contextPtr, maxSize );

    sendValues.reserve( contextPtr, maxSize );

    SCAI_LOG_DEBUG( logger, comm << ": recv buffer = " << recvValues << ", send buffer = " << sendValues );

    for ( PartitionId p = 0; p < nParts; ++p )
    {
        // compute the owner of the values that are in the current send buffer

        PartitionId actualPartition = comm.getNeighbor( -p );

        SCAI_LOG_INFO( logger, comm << ": will compute part for partition " << actualPartition << ", mData = " << *mData[actualPartition] )

        if ( p == 0 )
        {
            // start here with computation of own part

            SCAI_LOG_INFO( logger, comm << ": localX = " << localX << ", localY = " << localY )

            mData[actualPartition]->vectorTimesMatrix( sendValues, alphaValue, localX, betaValue, localY );
        }
        else
        {
            // This processor computes the part for actual partition and adds it

            SCAI_LOG_INFO( logger, comm << ": localX = " << localX << ", sendValues = " << sendValues )

            mData[actualPartition]->vectorTimesMatrix( sendValues, alphaValue, localX, ValueType( 1 ), sendValues );
        }

        SCAI_LOG_INFO( logger, comm << ": computed part for partition " << actualPartition << ", is " << sendValues );

        comm.shiftArray( recvValues, sendValues, COMM_DIRECTION );

        SCAI_LOG_DEBUG( logger, comm << ": received next " << recvValues );

        std::swap( sendValues, recvValues );
    }

    // we do not swap here as allocated data for localResult fits best

    utilskernel::HArrayUtils::assign( localResult, sendValues );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixPlusMatrix(
    const Scalar alpha,
    const _Matrix& matA,
    const Scalar beta,
    const _Matrix& matB )
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

    SCAI_LOG_INFO( logger, "Mat.plusMatrix, this = " << alpha << " * A + " << beta << " * B"
                   << ", A = " << A << ", B = " << B )

    // allocate this result matrix, but only if it is not aliased with A or B

    if ( this != &A && this != &B )
    {
        allocate( A.getRowDistributionPtr(), A.getColDistributionPtr() );
    }

    // Add matrices of each chunk

    SCAI_LOG_DEBUG( logger, "Mat.plusMatrix, mDataSize = " << mData.size() );

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->matrixPlusMatrix( alpha, *A.mData[i], beta, *B.mData[i] );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesMatrix(
    _Matrix& result,
    const Scalar alpha,
    const _Matrix& B,
    const Scalar beta,
    const _Matrix& C ) const
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
//We are calculating with a replicated _Matrix. So there is no need for an asyncronous call,
//because we have to sync in this method anyway (returning void not SyncToken)
// Note: any alias will be resolved by the matrix storage routine and not here
//       as it might introduce a temporary in any case
    res->mData[0]->matrixTimesMatrix( alpha.getValue<ValueType>(), *mData[0], *Bp->mData[0], beta.getValue<ValueType>(),
                                      *Cp->mData[0] );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseMatrix<ValueType>::maxNorm() const
{
    NormType<ValueType> myMaxDiff = 0;

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        NormType<ValueType> maxDiff = mData[i]->maxNorm();

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
NormType<ValueType> DenseMatrix<ValueType>::l1Norm() const
{
    const Communicator& comm = getRowDistribution().getCommunicator();

    NormType<ValueType> mySum = 0;
    IndexType n = mData.size();

    for ( IndexType i = 0; i < n; i++ )
    {
        mySum += static_cast<NormType<ValueType> >( mData[i]->l1Norm() );
    }

    return comm.sum( mySum );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseMatrix<ValueType>::l2Norm() const
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

    return common::Math::sqrt( comm.sum( mySum ) );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseMatrix<ValueType>::maxDiffNorm( const Matrix<ValueType>& other ) const
{
    if ( !( ( getNumColumns() == other.getNumColumns() ) && ( getNumRows() == other.getNumRows() ) ) )
    {
        COMMON_THROWEXCEPTION( "maxDiffNorm requires matrices of same format" );
    }

    // Implementation works only for same distributions and same type

    if ( ( getRowDistribution() == other.getRowDistribution() ) && ( getColDistribution() == other.getColDistribution() )
            && ( MatrixKind::DENSE == other.getMatrixKind() ) )
    {
        const DenseMatrix<ValueType>* typedOther = dynamic_cast<const DenseMatrix<ValueType>*>( &other );
        SCAI_ASSERT_DEBUG( typedOther, "SERIOUS: wrong dynamic cast: " << other )
        return maxDiffNormImpl( *typedOther );
    }
    else
    {
        SCAI_UNSUPPORTED( "maxDiffNorm requires temporary of " << other )
        DenseMatrix<ValueType> typedOther( other, getRowDistributionPtr(), getColDistributionPtr() );
        return maxDiffNormImpl( typedOther );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseMatrix<ValueType>::maxDiffNormImpl( const DenseMatrix<ValueType>& other ) const
{
    // implementation only supported for same distributions
    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution(), other.getRowDistribution() )
    SCAI_ASSERT_EQUAL_ERROR( getColDistribution(), other.getColDistribution() )

    typedef typename common::TypeTraits<ValueType>::AbsType AbsType;

    AbsType myMaxDiff = 0;

    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        AbsType maxDiff = mData[i]->maxDiffNorm( *other.mData[i] );

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

    // take the column data chunk that is owned by this processor regarding col dist

    const PartitionId myRank = getColDistribution().getCommunicator().getRank();
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

    // take the column data chunk that is owned by this processor regarding col dist

    const PartitionId myRank = getColDistribution().getCommunicator().getRank();
    return *mData[myRank];
}

template<typename ValueType>
IndexType DenseMatrix<ValueType>::getLocalNumValues() const
{
    // only locally stored number of values
    return getRowDistribution().getLocalSize() * getNumColumns();
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
    common::ScalarType type = common::getScalarType<ValueType>();
    stream << "DenseMatrix<" << type << ">( size = " << getNumRows() << " x " << getNumColumns() << ", rowdist = "
           << getRowDistribution() << ", coldist = " << getColDistribution() << " )";
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
    std::unique_ptr<DenseMatrix<ValueType> > newDenseMatrix( new DenseMatrix<ValueType>() );
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
_Matrix* DenseMatrix<ValueType>::create()
{
    return new DenseMatrix<ValueType>();
}

template<typename ValueType>
MatrixCreateKeyType DenseMatrix<ValueType>::createValue()
{
    common::ScalarType skind = common::getScalarType<ValueType>();
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

SCAI_COMMON_INST_CLASS( DenseMatrix, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */

