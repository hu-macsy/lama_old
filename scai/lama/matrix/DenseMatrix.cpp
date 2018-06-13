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
 * @author Thomas Brandes
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


// Metaprogramming to translate assign( _Matrix ) to assignImpl( Matrix<ValueType> )

template<typename ValueType, typename TList>
struct DenseMatrixWrapper;

template<typename ValueType>
struct DenseMatrixWrapper<ValueType, common::mepr::NullType>
{
    static void assign( DenseMatrix<ValueType>&, const _Matrix& other )
    {
        COMMON_THROWEXCEPTION( "DenseMatrix::assing: type of other matrix not supported --> " << other )
    }
};

template<typename ValueType, typename H, typename T>
struct DenseMatrixWrapper<ValueType, common::mepr::TypeList<H, T> >
{
    static void assign( DenseMatrix<ValueType>& obj, const _Matrix& other )
    {
        if ( other.getValueType() == common::getScalarType<H>() )
        {
            obj.assignImpl( static_cast<const Matrix<H>& >( other ) );
        }
        else
        {
            DenseMatrixWrapper<ValueType, T>::assign( obj, other );
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
DenseMatrix<ValueType>::DenseMatrix( ContextPtr ctx )
{
    computeOwners();
    allocateData( ctx );      // will initialize it with zero
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const IndexType numRows, const IndexType numColumns, ContextPtr ctx ) :

    Matrix<ValueType>( numRows, numColumns )

{
    computeOwners();
    allocateData( ctx );      // will initialize it with zero
}

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DistributionPtr rowDist, DistributionPtr colDist, ContextPtr ctx ) :

    Matrix<ValueType>( rowDist, colDist )

{
    computeOwners();
    allocateData( ctx );      // will initialize it with zero
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix(
    DistributionPtr rowDistribution,
    DenseStorage<ValueType> localStorage ) :

    Matrix<ValueType>( rowDistribution, std::make_shared<NoDistribution>( localStorage.getNumColumns() ) )
{
    // make some 'global' checks to verify correct sizes on all processors
    
    _Matrix::checkLocalStorageSizes( localStorage, *rowDistribution );

    mData.clear();
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( std::move( localStorage ) ) );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DenseStorage<ValueType> globalData ) :

    Matrix<ValueType>( globalData.getNumRows(), globalData.getNumColumns() )

{
    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( std::move( globalData ) ) );
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
DenseMatrix<ValueType>& DenseMatrix<ValueType>::operator=( DenseMatrix<ValueType>&& other )
{
    _Matrix::moveImpl( std::move( other ) );   // move sizes / distributions

    mData = std::move( other.mData );

    SCAI_ASSERT_EQ_ERROR( 0, other.mData.size(), "move of std::vector did not reset other to zero" )

    return *this;
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( const DenseMatrix<ValueType>& other ) :

    Matrix<ValueType>( other.getRowDistributionPtr(), other.getColDistributionPtr() )

{
    mData.resize ( other.mData.size() );

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        // due to copy we copy also the context from other

        mData[i].reset( other.mData[i]->copy() );
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
DenseMatrix<ValueType>::DenseMatrix( DenseMatrix<ValueType>&& other ) noexcept :

    Matrix<ValueType>( other ),
    mData( std::move( other.mData ) )

{
    // other might be inconsistent
    other._Matrix::setReplicatedMatrix( 0, 0 );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setIdentity( DistributionPtr dist )
{
    _Matrix::setDistributedMatrix( dist, dist );

    computeOwners();
    allocateData( getContextPtr() );   // initialized with zero

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
            ValueType one = 1;
            mData[i]->setDiagonal( one );
        }
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignDiagonal( const Vector<ValueType>& diagonal )
{
    DistributionPtr dist = diagonal.getDistributionPtr();

    _Matrix::setDistributedMatrix( dist, dist );

    computeOwners();  
    allocateData( getContextPtr() );   // initialized with zero

    const Communicator& comm = dist->getCommunicator();
    IndexType rank = comm.getRank();
    IndexType size = comm.getSize();
    SCAI_ASSERT_EQUAL_DEBUG( size, static_cast<IndexType>( mData.size() ) )

    if ( diagonal.getVectorKind() == VectorKind::DENSE )
    {
        const auto& denseDiagonal = static_cast<const DenseVector<ValueType>&>( diagonal );
        mData[rank]->assignDiagonal( denseDiagonal.getLocalValues() );
    }
    else
    {
        HArray<ValueType> localArray;
        diagonal.buildLocalValues( localArray );
        mData[rank]->assignDiagonal( localArray );
    }
}

/* ------------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::setCSRData(
    DistributionPtr rowDist,
    DistributionPtr colDist,
    const IndexType,
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
    mData[0]->setCSRData( n, m, ia, ja, values );

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
    mData[0]->setCSRData( rowIA.size() + 1, getNumColumns(), rowIA, rowJA, rowValues );
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
    allocateData( getContextPtr() );

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

    if ( other.getMatrixKind() == MatrixKind::DENSE && other.getValueType() == getValueType() )
    {
        assignTransposeImpl( static_cast<const DenseMatrix<ValueType>&>( other ) );
    }
    else
    {
        SCAI_UNSUPPORTED( "dense.assignTranspe( other ), other converted to DenseMatrix<" << getValueType() << ">" )
 
        assignTransposeImpl( convert<DenseMatrix<ValueType>>( other ) );
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::assignTransposeImpl( const DenseMatrix<ValueType>& matrix )
{
    SCAI_REGION( "Mat.Dense.transpose" )

    auto colDist = matrix.getColDistributionPtr();
    auto rowDist = matrix.getRowDistributionPtr();

    SCAI_LOG_INFO( logger, "transpose dense matrix with same value type, switch row/col distributions" )

    if ( rowDist->isReplicated() && colDist->isReplicated() )
    {
        DenseMatrix<ValueType> newMatrix( colDist, rowDist );

        newMatrix.mData[0]->assignTranspose( matrix.getLocalStorage() );

        *this = std::move( newMatrix );
    }
    else if ( rowDist->isReplicated() )
    {
        COMMON_THROWEXCEPTION( "transpose not supported for replicated matrices with distributed columns, matrix = " << matrix )
    }
    else if ( colDist->isReplicated() )
    {
        COMMON_THROWEXCEPTION( "transpose not supported for distributed matrices with replicated columns, matrix = " << matrix )
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( rowDist->getCommunicator(), colDist->getCommunicator(), "transpose only on same set of processors" )

        const Communicator& comm = rowDist->getCommunicator();

        const IndexType size = comm.getSize();

        DenseMatrix<ValueType> newMatrix( colDist, rowDist );

        // preparation for all2allv
 
        std::vector<IndexType> receiveSizes( size );
        std::vector<ValueType*> recvBuffer( size );

        for ( IndexType i = 0; i < size; ++i )
        {
            HArray<ValueType>& recvData = newMatrix.mData[i]->getData();

            recvBuffer[i]   = hmemo::hostWriteAccess( recvData ).get();
            receiveSizes[i] = recvData.size();
        }

        std::vector<IndexType> sendSizes( size );
        std::vector<const ValueType*> sendBuffer( size );

        for ( IndexType i = 0; i < size; ++i )
        {
            const HArray<ValueType>& sendData = matrix.mData[i]->getValues();

            sendSizes[i] = sendData.size();
            sendBuffer[i] = hmemo::hostReadAccess( sendData ).get();
        }

        // MPI call

        comm.all2allv( recvBuffer.data(), receiveSizes.data(), sendBuffer.data(), sendSizes.data() );

        for ( IndexType i = 0; i < size; ++i )
        {
            DenseStorage<ValueType>& storage = *newMatrix.mData[i];
            HArray<ValueType>& data = storage.getData();
            utilskernel::HArrayUtils::transpose( data, storage.getNumRows(), storage.getNumColumns(), data, false );
        }

        *this = std::move( newMatrix );
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assign( const _Matrix& other )
{
    SCAI_LOG_INFO( logger, "assign " << other << " to " << *this )

    if ( &other == this )
    {
        SCAI_LOG_INFO( logger, "self assign, is skpped" )
    }
    else 
    {
        DenseMatrixWrapper<ValueType, SCAI_NUMERIC_TYPES_HOST_LIST>::assign( *this, other );
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::assignImpl( const Matrix<OtherValueType>& other )
{
    SCAI_ASSERT_ERROR( getContextPtr(), "assign to matrix without context" )

    if ( other.getMatrixKind() == MatrixKind::DENSE )
    {
        assignDense( static_cast<const DenseMatrix<OtherValueType>&>( other ) );
    }
    else if ( other.getMatrixKind() == MatrixKind::SPARSE )
    {
        assignSparse( static_cast<const SparseMatrix<OtherValueType>&>( other ) );
    }
    else
    {
        SCAI_THROWEXCEPTION( common::InvalidArgumentException, "Unsupported matrix kind, other = " << other )
    }
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::assignSparse( const SparseMatrix<OtherValueType>& other )
{
    SCAI_LOG_INFO( logger, "assignSparse: other = " << other << " @ " << getContextPtr() )

    // we need replicated column distribution to get this routine working

    if ( !other.getColDistribution().isReplicated() )
    {
        auto repColDist = std::make_shared<NoDistribution>( other.getNumColumns() );
        auto repOther   = distribute<CSRSparseMatrix<ValueType>>( other, other.getRowDistributionPtr(), repColDist );

        assignSparse( repOther );

        splitColumns( other.getColDistributionPtr() );

        return;
    }

    // replicated columns in sparse matrix, so we can assign local data

    _Matrix::setDistributedMatrix( other.getRowDistributionPtr(), other.getColDistributionPtr() );

    ContextPtr ctx = getContextPtr();

    mData.resize( 1 );
    mData[0].reset( new DenseStorage<ValueType>( ctx ) );
    mData[0]->assign( other.getLocalStorage() );

    computeOwners();
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void DenseMatrix<ValueType>::assignDense( const DenseMatrix<OtherValueType>& other )
{
    // check for valid pointer, might be dynamic cast went wrong somewhere else
    //SCAI_ASSERT_ERROR( &other, "NULL matrix in assignment operator" )
    SCAI_LOG_INFO( logger, "assign dense, this = " << this << ", other = " << &other )
    // inherit size and distributions
    _Matrix::setDistributedMatrix( other.getRowDistributionPtr(), other.getColDistributionPtr() );
    mData.resize( other.mData.size() );
    IndexType n = static_cast<IndexType>( other.mData.size() );

    for ( IndexType i = 0; i < n; ++i )
    {
        SCAI_LOG_DEBUG( logger, "copy block " << i << " of " << n << " = " << *other.mData[i] )
        mData[i].reset( new DenseStorage<ValueType>() );
        mData[i]->assign( *other.mData[i] );
    }
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
    mData[0].reset( new DenseStorage<ValueType>( getContextPtr() ) );
    mData[0]->assign( storage );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignLocal( const _MatrixStorage& storage, DistributionPtr rowDist )
{
    _Matrix::checkLocalStorageSizes( storage, *rowDist );

    IndexType numColumns = storage.getNumColumns();  // same for all processors

    _Matrix::setDistributedMatrix( rowDist, std::make_shared<NoDistribution>( numColumns ) );

    mData.resize( 1 );

    if ( !mData[0] )
    {
        mData[0].reset( new DenseStorage<ValueType>() );
    }

    mData[0]->assign( storage );
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::assignDistribute( const _MatrixStorage& storage, DistributionPtr rowDist, DistributionPtr colDist )
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

            splitColumnData( mData, static_cast<const DenseStorage<ValueType>&>( storage ), *colDist );
        }
        else if ( colDist->isReplicated() )
        {
            mData.resize( 1 );
            mData[0].reset( new DenseStorage<ValueType>() );
            mData[0]->assign( storage );
        }
        else
        {
            // conversion to dense storage needed before splitting

            auto localData = convert<DenseStorage<ValueType>>( storage );
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
void DenseMatrix<ValueType>::assignDistribute( const _Matrix& other, DistributionPtr rowDist, DistributionPtr colDist )
{
    assign( other );
    redistribute( rowDist, colDist );
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

        HArray<ValueType> denseData;   // output array for joined column data

        joinColumnData( denseData, 0, numLocalRows );

        DenseStorage<ValueType> denseStorage( numLocalRows, getNumColumns(), std::move( denseData ) );

        storage = denseStorage;
    }

    SCAI_LOG_DEBUG( logger, "buildLocalStorage( " << *this << " ) = " << storage )
}

/* ------------------------------------------------------------------ */

template<typename ValueType>
void DenseMatrix<ValueType>::disassemble( 
    MatrixAssembly<ValueType>& assembly,
    const IndexType rowOffset,
    const IndexType colOffset ) const
{
    SCAI_LOG_INFO( logger, "build local storage with replicated columns for " << *this )

    const IndexType numLocalRows = getRowDistribution().getLocalSize();

    HArray<IndexType> rowLocal2Global; 
    HArray<IndexType> offsets;
    HArray<IndexType> colLocal2Global;
  
    getRowDistribution().getOwnedIndexes( rowLocal2Global );
    getColDistribution().getAnyLocal2Global( offsets, colLocal2Global );

    auto rOffsets         = hostReadAccess( offsets );
    auto rColLocal2Global = hostReadAccess( colLocal2Global );
    auto rRowLocal2Global = hostReadAccess( rowLocal2Global );

    ValueType zero = 0;

    for ( size_t p = 0; p < mData.size(); ++p )
    {
        auto values = hostReadAccess( mData[p]->getValues() );

        IndexType offset          = rOffsets[p];
        IndexType numLocalColumns = rOffsets[p + 1] - rOffsets[p];

        for ( IndexType localJ = 0; localJ < numLocalColumns; ++localJ )
        {
            IndexType globalJ = rColLocal2Global[ offset + localJ ];

            for ( IndexType i = 0; i < numLocalRows ; ++i )
            {
                ValueType val = values[ localJ + i * numLocalColumns ];
                
                if ( val == zero )
                {
                    continue;
                }

                IndexType globalI = rRowLocal2Global[i];

                assembly.push( globalI + rowOffset, globalJ + colOffset, val );
            }
        }
    }
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
            ReadAccess<ValueType> rLocalData( mData[p]->getValues(), hostContext );

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
void DenseMatrix<ValueType>::allocateData( ContextPtr ctx )
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
            mData[0].reset( new DenseStorage<ValueType>( numLocalRows, getNumColumns(), ctx ) );
        }

        return;
    }

    IndexType count = 0;   // sum up the sizes, verify correct sum

    for ( PartitionId p = 0; p < numChunks; ++p )
    {
        IndexType numLocalColumns = colDist.getAnyLocalSize( p );
        count += numLocalColumns;
        mData[p].reset( new DenseStorage<ValueType>( numLocalRows, numLocalColumns, ctx ) );
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

    ReadAccess<ValueType> columnDataRead( columnData.getValues(), ctx );
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
        HArray<ValueType> colData;  // output array for joined column data, moved into dense storage
        joinColumnData( colData, 0, numLocalRows );
        mData.clear();
        mData.resize( 1 );
        mData[0].reset( new DenseStorage<ValueType>( numLocalRows, numCols, std::move( colData ) ) );
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
void DenseMatrix<ValueType>::resize( DistributionPtr rowDistributionPtr, DistributionPtr colDistributionPtr )
{
    // disassemble this matrix and fill it up again

    MatrixAssembly<ValueType> assembly;

    this->disassemble( assembly );

    IndexType newNumRows = rowDistributionPtr->getGlobalSize();
    IndexType newNumCols = colDistributionPtr->getGlobalSize();

    // truncate elements if necessary to avoid ERROR messages when we fil up

    if ( newNumRows < getNumRows() || newNumCols < getNumColumns() )
    {
        assembly.truncate( newNumRows, newNumCols );
    }

    // and now fill the assembly back

    allocate( rowDistributionPtr, colDistributionPtr );

    this->fillFromAssembly( assembly );
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
    ReadAccess<ValueType> repData( global.getValues(), contextPtr );
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
    ReadAccess<ValueType> localVals( distributedData.getValues(), contextPtr );
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
    redistributor.redistributeN( newLocalData.getData(), oldLocalData.getValues(), nCols );
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
        getLocalStorage().getRow( row, localRowIndex );
        return;
    }

    // with column distribution: join the corresponding column data

    joinColumnData( row, localRowIndex, 1 );

    SCAI_LOG_INFO( logger, "local row with joined column data: " << row )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getRowLocal( Vector<ValueType>&, const IndexType ) const
{
    COMMON_THROWEXCEPTION( "not available yet" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getRow( Vector<ValueType>& row, const IndexType globalRowIndex ) const
{
    // if v is not a dense vector or not of same type, use a temporary dense vector

    if ( row.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_LOG_WARN( logger, "getRow requires temporary" )
        DenseVector<ValueType> denseRow;
        getRow( denseRow, globalRowIndex );
        row.assign( denseRow );   // transform the dense vector into sparse vector
        return;
    }

    SCAI_REGION( "Mat.Dense.getRow" )

    DenseVector<ValueType>& denseRow = static_cast<DenseVector<ValueType>&>( row );

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
        auto recvPlan = CommunicationPlan::buildBySizes( NULL, 0 );  // nothing to receive

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

        auto sendPlan = CommunicationPlan::buildBySizes( NULL, 0 );
        auto recvPlan = CommunicationPlan::buildBySizes( NULL, 0 );

        recvPlan.singleEntry( rowOwner, size );

        SCAI_LOG_DEBUG( logger, comm << ": getRow, recvPlan = " << recvPlan << ", sendPlan = " << sendPlan )

        comm.exchangeByPlan( values, recvPlan, dummySend, sendPlan );
    }

    // guarantee consistency in the dense vector for the local data

    SCAI_ASSERT_EQ_ERROR( values.size(), getColDistribution().getLocalSize(), "serious mismatch" )
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getColumn( Vector<ValueType>& col, const IndexType globalColIndex ) const
{
    // if col is not a dense vector, use a temporary dense vector

    if ( col.getVectorKind() != VectorKind::DENSE )
    {
        SCAI_LOG_WARN( logger, "getCol requires temporary, use DenseVector on DenseMatrix" )
        DenseVector<ValueType> denseColumn;
        getColumn( denseColumn, globalColIndex );
        col.assign( denseColumn );   // transform the dense vector into sparse vector, works for all
        return;
    }

    SCAI_REGION( "Mat.Dense.getColumn" )

    SCAI_ASSERT_DEBUG( dynamic_cast<DenseVector<ValueType>*>( &col ), "col not DenseVector<" << getValueType() << ">" )

    DenseVector<ValueType>& denseCol = static_cast<DenseVector<ValueType>&>( col );

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

        mData[0]->setRow( row, localRowIndex, op );
        return;
    }

    HArray<IndexType> offsets;
    HArray<IndexType> perm;

    getColDistribution().getAnyLocal2Global( offsets, perm );

    HArray<ValueType> rowResorted;   // row resorted according to the owners

    utilskernel::HArrayUtils::gather( rowResorted, row, perm, common::BinaryOp::COPY );

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

        mData[ip]->setRow( rowPartition, localRowIndex, op );
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

    mData[owner]->setColumn( column, localColIndex, op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::getDiagonal( Vector<ValueType>& diagonal ) const
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( diagonal.getVectorKind() != VectorKind::DENSE )
    {
        // MIGHT BE WORTH A WARNING
        DenseVector<ValueType> denseDiagonal;
        getDiagonal( denseDiagonal );
        diagonal = denseDiagonal;
        return;
    }

    // we can recast it now to dense vector, so we have access to its local values

    DenseVector<ValueType>& diagonalDense = static_cast<DenseVector<ValueType>&>( diagonal );
    diagonalDense.allocate( getRowDistributionPtr() );
    getLocalStorage().getDiagonal( diagonalDense.getLocalValues() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const Vector<ValueType>& diagonal )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "setDiagonal only for square matrices with same row/col distribution" )
    }

    if ( getRowDistribution() != diagonal.getDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    if ( diagonal.getVectorKind() != VectorKind::DENSE )
    {
        // MIGHT BE WORTH A WARNING
        setDiagonal( convert<DenseVector<ValueType>>( diagonal ) );
        return;
    }

    const DenseVector<ValueType>& diagonalDense = static_cast<const DenseVector<ValueType>&>( diagonal );

    getLocalStorage().setDiagonalV( diagonalDense.getLocalValues() );
}

template<typename ValueType>
void DenseMatrix<ValueType>::setDiagonal( const ValueType& diagonalValue )
{
    if ( getRowDistribution() != getColDistribution() )
    {
        COMMON_THROWEXCEPTION( "Diagonal calculation only for equal distributions." )
    }

    getLocalStorage().setDiagonal( diagonalValue );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::reduce(
    Vector<ValueType>& v, 
    const IndexType dim, 
    const common::BinaryOp reduceOp, 
    const common::UnaryOp elemOp ) const
{
    if ( v.getVectorKind() == VectorKind::DENSE )
    {
        reduceImpl( static_cast<DenseVector<ValueType>&>( v ), dim, reduceOp, elemOp );
    }
    else
    {
        DenseVector<ValueType> tmpV;
        reduceImpl( tmpV, dim, reduceOp, elemOp );
        v = tmpV;
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::reduceImpl(
    DenseVector<ValueType>& v, 
    const IndexType dim, 
    const common::BinaryOp reduceOp, 
    const common::UnaryOp elemOp ) const
{
    SCAI_REGION( "Mat.Dense.reduce" )

    if ( dim == 0 )
    {
        v.allocate( getRowDistributionPtr() );

        v = ValueType( 0 );   // initialize v with neutral element

        for ( size_t k = 0; k < mData.size(); ++k )
        {
            mData[k]->reduce( v.getLocalValues(), 0, reduceOp, elemOp );
        }

        return;
    }

    if ( dim == 1 )
    {
        v.allocate( getColDistributionPtr() );

        v = ValueType( 0 );   // initialize v with neutral element

        if ( getRowDistribution().getCommunicator().getSize() == 1 )
        {
            // full matrix is replicated, the columns might have any distribution

            PartitionId rank = getColDistribution().getCommunicator().getRank();

            mData[rank]->reduce( v.getLocalValues(), 1, reduceOp, elemOp );

            return;   // matrix is replicated, compute just my values
        }

        // rows are distributed

        IndexType np = getColDistribution().getCommunicator().getSize();

        if ( np == 1 )
        {
             SCAI_ASSERT_EQ_ERROR( reduceOp, common::BinaryOp::ADD, "only add supported" )

             mData[0]->reduce( v.getLocalValues(), 1, reduceOp, elemOp );
             getRowDistribution().getCommunicator().sumArray( v.getLocalValues() );
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

        utilskernel::HArrayUtils::assign( sendValues, v.getLocalValues() );

        const Communicator& comm = getColDistribution().getCommunicator();

        for ( PartitionId p = 0; p < np; ++p )
        {
            // compute the owner of the values that are in the current send buffer
            PartitionId actualPartition = comm.getNeighbor( -p );
            mData[actualPartition]->reduce( sendValues, 1, reduceOp, elemOp );
            comm.shiftArray( recvValues, sendValues, COMM_DIRECTION );
            std::swap( sendValues, recvValues );
        }

        utilskernel::HArrayUtils::assign( v.getLocalValues(), sendValues );
    }
    else
    {
        COMMON_THROWEXCEPTION( "illegal reduce dim = " << dim << " for dense matrix" )
    }
}

/* -------------------------------------------------------------------------- */
/*   scaling of matrix entries                                                */
/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::scaleRows( const DenseVector<ValueType>& scaleY )
{
    SCAI_ASSERT_EQ_ERROR( getRowDistribution(), scaleY.getDistribution(), 
                          "distribution of scale vector does not match" )

    const HArray<ValueType>& localY = scaleY.getLocalValues();

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->scaleRows( localY );
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::scale( const ValueType& alpha )
{
    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->scale( alpha );
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::conj()
{
    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->conj();
    }
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
        IndexType  jLocal = invalidIndex;

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

        SCAI_ASSERT_ERROR( jLocal != invalidIndex, "non local column index" )
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

    if ( iLocal == invalidIndex )
    {
        return; // this processor does not have the value
    }

    const Distribution& distributionCol = getColDistribution();

    PartitionId owner  = 0;

    IndexType   jLocal = invalidIndex;

    if ( distributionCol.getNumPartitions() == 1 )
    {
        jLocal = j;
    }
    else 
    {
        owner  = distributionCol.getAnyOwner( j );
        jLocal = distributionCol.getAnyLocalIndex( j, owner );
    }

    SCAI_ASSERT_ERROR( jLocal != invalidIndex, "non local column index" )

    mData[owner]->setValue( iLocal, jLocal, val, op );
}

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesScalar( const Matrix<ValueType>& other, ValueType alpha )
{
    SCAI_LOG_INFO( logger, " this = " << alpha << " * " << other )
    assign( other );
    SCAI_LOG_INFO( logger, " this = other = " << *this )

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->scale( alpha );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesVectorDense(
    DenseVector<ValueType>& result,
    const ValueType alpha,
    const DenseVector<ValueType>& x,
    const ValueType beta,
    const DenseVector<ValueType>* y,
    const common::MatrixOp op ) const
{
    if ( common::MatrixOp::NORMAL == op )
    {
        matrixTimesVectorImpl( result, alpha, x, beta, y );
    }
    else if ( common::MatrixOp::TRANSPOSE != op )
    {
        COMMON_THROWEXCEPTION( op << " not supported yet for matrix * vector" )
    }
    else if ( this->getColDistribution().getCommunicator().getSize() == 1 )
    {
        // Each processor has full columns, resultVector is replicated, communication only needed to sum up results
        // use routine provided by this CRTP

        Matrix<ValueType>::vectorTimesMatrixRepCols( result, alpha, x, beta, y );
    }
    else
    {
        vectorTimesMatrixImpl( result, alpha, x, beta, y );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesVectorImpl(
    DenseVector<ValueType>& denseResult,
    const ValueType alphaValue,
    const DenseVector<ValueType>& denseX,
    const ValueType betaValue,
    const DenseVector<ValueType>* denseY ) const
{
    SCAI_REGION( "Mat.Dense.timesVector" )

    HArray<ValueType>& localResult = denseResult.getLocalValues();
    const HArray<ValueType>& localY = denseY == nullptr ? localResult : denseY->getLocalValues();
    ContextPtr localContext = mData[0]->getContextPtr();
    const Distribution& colDist = getColDistribution();
    const Communicator& comm = colDist.getCommunicator();
    PartitionId rank = comm.getRank();
    PartitionId n = colDist.getNumPartitions();
    mData[0]->prefetch();

    // It makes no sense to prefetch denseX because, if a transfer is started
    // the halo update needs to wait for this transfer to finish

    if ( denseY != nullptr )
    {
        denseY->prefetch( localContext );
    }

    const HArray<ValueType>& localX = denseX.getLocalValues();

    SCAI_LOG_INFO( logger,
                   comm << ": matrixTimesVector" << ", alpha = " << alphaValue << ", localX = " << localX << ", beta = " << betaValue << ", localY = " << localY )
    SCAI_LOG_INFO( logger,
                   "Aliasing: result = y : " << ( &denseResult == denseY ) << ", local = " << ( &localResult == &localY ) )

    if ( n == 1 )
    {
// replicated column distribution, only on local block, X is replicated
// localResult = alpha * mData[0] * X + beta * localY
        const DenseStorage<ValueType>& dense = *mData[0];
        SCAI_LOG_INFO( logger, comm << ": matrixTimesVector, single dense block = " << dense )
        dense.matrixTimesVector( localResult, alphaValue, localX, betaValue, localY, common::MatrixOp::NORMAL );
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

    if ( SyncKind::SYNCHRONOUS != _Matrix::getCommunicationKind() )
    {
        SCAI_LOG_INFO( logger, comm << ": asynchronous communication" )
// asynchronous communication always requires same sizes of arrays, might shift some more data
        std::unique_ptr<tasking::SyncToken> st( comm.shiftAsync( *recvValues, *sendValues, COMM_DIRECTION ) );
        SCAI_LOG_INFO( logger,
                       comm << ": matrixTimesVector, my dense block = " << *mData[rank] << ", localX = " << localX << ", localY = " << localY << ", localResult = " << localResult )
// overlap communication with local computation
        mData[rank]->matrixTimesVector( localResult, alphaValue, localX, betaValue, localY, common::MatrixOp::NORMAL );
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
            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, x, static_cast<ValueType>( 1.0 ), localResult, common::MatrixOp::NORMAL );
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
        mData[rank]->matrixTimesVector( localResult, alphaValue, localX, betaValue, localY, common::MatrixOp::NORMAL );
        std::swap( sendValues, recvValues );

        for ( PartitionId p = 1; p < n; ++p )
        {
            PartitionId actualPartition = comm.getNeighbor( -p );
            comm.shiftArray( *recvValues, *sendValues, COMM_DIRECTION );
            SCAI_LOG_DEBUG( logger,
                            comm << ": send " << *sendValues << ", recv " << *recvValues << ", actual = " << actualPartition )
            SCAI_LOG_INFO( logger,
                           comm << ": matrixTimesVector, actual dense block [" << actualPartition << "] = " << *mData[actualPartition] << ", sendX = " << *sendValues << ", localResult = " << localResult )
            mData[actualPartition]->matrixTimesVector( localResult, alphaValue, *sendValues, ValueType( 1 ), localResult, common::MatrixOp::NORMAL );
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
    const DenseVector<ValueType>* denseY ) const
{
    SCAI_REGION( "Mat.Dense.vectorTimesMatrix" )

    HArray<ValueType>& localResult = denseResult.getLocalValues();

    const HArray<ValueType>& localY = denseY == nullptr ? localResult : denseY->getLocalValues();

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
        mData[rank]->matrixTimesVector( localResult, alphaValue, localX, betaValue, localY, common::MatrixOp::TRANSPOSE );

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

            mData[actualPartition]->matrixTimesVector( sendValues, alphaValue, localX, betaValue, localY, common::MatrixOp::TRANSPOSE );
        }
        else
        {
            // This processor computes the part for actual partition and adds it

            SCAI_LOG_INFO( logger, comm << ": localX = " << localX << ", sendValues = " << sendValues )

            mData[actualPartition]->matrixTimesVector( sendValues, alphaValue, localX, ValueType( 1 ), sendValues, common::MatrixOp::TRANSPOSE );
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
void DenseMatrix<ValueType>::binaryOp(
    const Matrix<ValueType>& matrixA,
    const common::BinaryOp op,
    const Matrix<ValueType>& matrixB )
{
    SCAI_LOG_INFO( logger, "this = " << "A " << op << " B, A = " << matrixA << ", B = " << matrixB )

    if ( matrixA.getMatrixKind() != MatrixKind::DENSE )
    {
        if ( &matrixB == this )
        {   
            auto denseA = convert<DenseMatrix<ValueType>>( matrixA );
            binaryOp( denseA, op, matrixB );
        }
        else
        {
            // reuse this storage for conversion of a
            assign( matrixA );
            binaryOp( *this, op, matrixB );
        }
    }
    else if ( matrixB.getMatrixKind() != MatrixKind::DENSE )
    {
        if ( &matrixA == this )
        {
            auto denseB = convert<DenseMatrix<ValueType>>( matrixB );
            binaryOp( matrixA, op, denseB );
        }
        else
        {
            // reuse this matrix for conversion of a
            assign( matrixB ); 
            binaryOp( matrixA, op, *this );
        }
    }
    else
    {
        // Here we can call binary op for dense matrices

        binaryOpDense( static_cast<const DenseMatrix<ValueType>&>( matrixA ), op,
                       static_cast<const DenseMatrix<ValueType>&>( matrixB ) );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const Matrix<ValueType>& matA,
    const ValueType beta,
    const Matrix<ValueType>& matB )
{
    SCAI_ASSERT_EQ_ERROR( matA.getRowDistribution(), matB.getRowDistribution(), "size/dist mismatch of matrices to add" )
    SCAI_ASSERT_EQ_ERROR( matB.getColDistribution(), matB.getColDistribution(), "size/dist mismatch of matrices to add" )

    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", A = " << matA << ", B = " << matB )

    SCAI_ASSERT_EQ_ERROR( matA.getMatrixKind(), MatrixKind::DENSE, "denseMatrix = alpha * matA + beta * matB, matA must be dense" )
    SCAI_ASSERT_EQ_ERROR( matB.getMatrixKind(), MatrixKind::DENSE, "denseMatrix = alpha * matA + beta * matB, matB must be dense" )

    const DenseMatrix<ValueType>& denseA = static_cast<const DenseMatrix<ValueType>&>( matA );
    const DenseMatrix<ValueType>& denseB = static_cast<const DenseMatrix<ValueType>&>( matB );

    // Now we can add dense matrices

    matrixPlusMatrixDense( alpha, denseA, beta, denseB );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixPlusMatrixDense(
    const ValueType alpha,
    const DenseMatrix<ValueType>& A,
    const ValueType beta,
    const DenseMatrix<ValueType>& B )
{
    SCAI_REGION( "Mat.plusMatrix" )

    // already verified

    SCAI_ASSERT_EQ_DEBUG( A.getRowDistribution(), B.getRowDistribution(), "size/dist mismatch of matrices to add" )
    SCAI_ASSERT_EQ_DEBUG( B.getColDistribution(), B.getColDistribution(), "size/dist mismatch of matrices to add" )

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
void DenseMatrix<ValueType>::binaryOpDense(
    const DenseMatrix<ValueType>& matrixA,
    const common::BinaryOp op,
    const DenseMatrix<ValueType>& matrixB )
{
    SCAI_REGION( "Mat.Dense.binaryOp" )

    SCAI_ASSERT_EQ_ERROR( matrixA.getRowDistribution(), matrixB.getRowDistribution(), "size/dist mismatch of matrices in binaryOp" )
    SCAI_ASSERT_EQ_ERROR( matrixB.getColDistribution(), matrixB.getColDistribution(), "size/dist mismatch of matrices in binaryOp" )

    // Now we can do it completely local

    if ( this != &matrixA && this != &matrixB )
    {
        allocate( matrixA.getRowDistributionPtr(), matrixA.getColDistributionPtr() );
    }

    // Add matrices of each chunk

    SCAI_LOG_DEBUG( logger, "Mat.plusMatrix, mDataSize = " << mData.size() );

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        mData[i]->binaryOp( *matrixA.mData[i], op, *matrixB.mData[i] );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::matrixTimesMatrix(
    Matrix<ValueType>& result,
    const ValueType alpha,
    const Matrix<ValueType>& B,
    const ValueType beta,
    const Matrix<ValueType>& C ) const
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
    else if ( res == Cp && beta != common::Constants::ZERO )
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
    res->mData[0]->matrixTimesMatrix( alpha, *mData[0], *Bp->mData[0], beta, *Cp->mData[0] );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseMatrix<ValueType>::selectComplexPart( Matrix<RealType<ValueType> >& x, common::ComplexPart kind ) const
{
    if ( kind == common::ComplexPart::REAL )
    {
        x = cast<RealType<ValueType>>( *this );
    }
    else
    {
        ValueType i = common::TypeTraits<ValueType>::imaginaryUnit();
        DenseMatrix<ValueType> tmp( *this );
        tmp *= -i;  // imaginary part becomes real part
        x = cast<RealType<ValueType>>( tmp );
    }
}

template<typename ValueType>
void DenseMatrix<ValueType>::buildComplex( const Matrix<RealType<ValueType> >& x, const Matrix<RealType<ValueType> >& y )
{
    SCAI_LOG_INFO( logger, "buildComplex<" << getValueType() << ">( x, y ) with x = " << x << ", y = " << y )

    auto x1 = convert<DenseMatrix<ValueType>>( x );
    auto y1 = convert<DenseMatrix<ValueType>>( y );
    ValueType i = common::TypeTraits<ValueType>::imaginaryUnit(); 
    matrixPlusMatrix( 1, x1, i, y1 );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseMatrix<ValueType>::maxNorm() const
{
    RealType<ValueType> myMaxDiff = 0;

    for ( size_t i = 0; i < mData.size(); ++i )
    {
        RealType<ValueType> maxDiff = mData[i]->maxNorm();

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
RealType<ValueType> DenseMatrix<ValueType>::l1Norm() const
{
    const Communicator& comm = getRowDistribution().getCommunicator();

    RealType<ValueType> mySum = 0;
    IndexType n = mData.size();

    for ( IndexType i = 0; i < n; i++ )
    {
        mySum += static_cast<RealType<ValueType> >( mData[i]->l1Norm() );
    }

    return comm.sum( mySum );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DenseMatrix<ValueType>::l2Norm() const
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
RealType<ValueType> DenseMatrix<ValueType>::maxDiffNorm( const Matrix<ValueType>& other ) const
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
        return maxDiffNormImpl( distribute<DenseMatrix<ValueType>>( other, getRowDistributionPtr(), getColDistributionPtr() ) );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseMatrix<ValueType>::maxDiffNormImpl( const DenseMatrix<ValueType>& other ) const
{
    // implementation only supported for same distributions
    SCAI_ASSERT_EQUAL_ERROR( getRowDistribution(), other.getRowDistribution() )
    SCAI_ASSERT_EQUAL_ERROR( getColDistribution(), other.getColDistribution() )

    typedef typename common::TypeTraits<ValueType>::RealType RealType;

    RealType myMaxDiff = 0;

    for ( unsigned int i = 0; i < mData.size(); ++i )
    {
        RealType maxDiff = mData[i]->maxDiffNorm( *other.mData[i] );

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
    newDenseMatrix->allocate( getRowDistributionPtr(), getColDistributionPtr() );
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

