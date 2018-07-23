/**
 * @file COOStorage.cpp
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
 * @brief Instantitions for template class COOStorage.
 * @author Lauretta Schubert, Thomas Brandes
 * @date 25.05.2011
 */

// hpp
#include <scai/lama/storage/COOStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/COOUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/hmemo.hpp>

#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/loop.hpp>
#include <scai/common/macros/unsupported.hpp>

// sqrt for all value types
#include <cmath>
#include <memory>

using std::unique_ptr;

using namespace scai::hmemo;

namespace scai
{

using tasking::SyncToken;

using utilskernel::HArrayUtils;

using sparsekernel::COOUtils;
using sparsekernel::CSRUtils;

using common::BinaryOp;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, COOStorage<ValueType>::logger, "MatrixStorage.COOStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( ContextPtr ctx ) : 

    MatrixStorage<ValueType>( 0, 0, ctx ),
    mIA( ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "COOStorage for matrix " << getNumRows() 
                             << " x " << getNumColumns() << ", no non-zero elements @ " << *ctx )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( IndexType numRows, IndexType numColumns, ContextPtr ctx ) : 

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mIA( ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "COOStorage for matrix " << getNumRows() 
                             << " x " << getNumColumns() << ", no non-zero elements @ " << *ctx )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::verifySorting()
{
    COOUtils::normalize( mIA, mJA, mValues, common::BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( 
    const IndexType numRows, 
    const IndexType numColumns,
    HArray<IndexType> ia, 
    HArray<IndexType> ja, 
    HArray<ValueType> values, 
    ContextPtr ctx ) : 

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),

    mIA( std::move( ia ) ),
    mJA( std::move( ja ) ),
    mValues( std::move( values ) )

{
    SCAI_ASSERT_EQ_ERROR( mIA.size(), mJA.size(), "serious mismatch for sizes of ia, ja" )
    SCAI_ASSERT_EQ_ERROR( mIA.size(), mValues.size(), "serious mismatch for sizes of ia, ja" )

    verifySorting();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( const COOStorage<ValueType>& other ) : 

    MatrixStorage<ValueType>( other )

{
    assignImpl( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( COOStorage<ValueType>&& other ) :

    MatrixStorage<ValueType>( std::move( other ) ),

    mIA( std::move( other.mIA ) ),
    mJA( std::move( other.mJA ) ),
    mValues( std::move( other.mValues ) )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::~COOStorage()
{
    SCAI_LOG_DEBUG( logger, "~COOStorage for matrix " << getNumRows() << " x " << getNumColumns() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>& COOStorage<ValueType>::operator=( const COOStorage<ValueType>& other )
{
    assignCOO( other );
    return *this;
}

template<typename ValueType>
COOStorage<ValueType>& COOStorage<ValueType>::operator=( COOStorage<ValueType>&& other )
{
    // call of move assignment for base class 

    MatrixStorage<ValueType>::moveImpl( std::move( other ) );

    // move of all member variables

    mIA = std::move( other.mIA );
    mJA = std::move( other.mJA );
    mValues = std::move( other.mValues );

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::StorageWrapper<COOStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    if ( other.getFormat() == Format::COO )
    {
        // both storage have COO format, use special method for it

        assignCOO( static_cast<const COOStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRDataImpl( otherCSR.getNumRows(), otherCSR.getNumColumns(), 
                        otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues() );
    }
    else 
    {
        HArray<IndexType>  csrIA;
        HArray<IndexType>  csrJA;
        HArray<ValueType>  csrValues;     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        // just a thought for optimization: use mIA, mJA, mValues instead of csrIA, csrJA, csrValues
        // but does not help much at all as resort of entries requires already temporaries.

        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::assignCOO( const COOStorage<OtherValueType>& other )
{
    if ( other.getValueType() == getValueType() )
    {
        if ( reinterpret_cast<const COOStorage<ValueType>*>( &other ) == this )
        {
            SCAI_LOG_INFO( logger, typeName() << ": self assign, skipped, matrix = " << other )
            return;
        }
    }

    auto ctx = getContextPtr();

    // both storage have COO format, we can just copy the corresponding arrays to the right context

    _MatrixStorage::_assign( other );

    HArrayUtils::setArray( mIA, other.getIA(), common::BinaryOp::COPY, ctx );
    HArrayUtils::setArray( mJA, other.getJA(), common::BinaryOp::COPY, ctx );
    HArrayUtils::setArray( mValues, other.getValues(), common::BinaryOp::COPY, ctx );

    // normalize not required as other storage should contain normalized data
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format COOStorage<ValueType>::getFormat() const
{
    return Format::COO;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "COOStorage " << getNumRows() << " x " << getNumColumns() << ", #values = " << getNumValues() << endl;
    ContextPtr host = Context::getHostPtr();
    ReadAccess<IndexType> ia( mIA, host );
    ReadAccess<IndexType> ja( mJA, host );
    ReadAccess<ValueType> values( mValues, host );

    for ( IndexType i = 0; i < mIA.size(); i++ )
    {
        stream << "@[ " << ia[i] << ", " << ja[i] << " ] = " << values[i] << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::clear()
{
    _MatrixStorage::setDimension( 0, 0 );

    mIA.clear();
    mJA.clear();
    mValues.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::check( const char* msg ) const
{
    SCAI_ASSERT_EQ_ERROR( mIA.size(), mJA.size(), msg << ": serious mismatch for sizes of ia, ja" )
    SCAI_ASSERT_EQ_ERROR( mIA.size(), mValues.size(), msg << ": serious mismatch for sizes of ia, values" )

    const IndexType m = getNumRows();
    const IndexType n = getNumColumns();

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mIA, m ), msg << ": invalid row indexes, #rows = " << m );
    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mJA, n ), msg << ": invalid column indexes, #cols = " << n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setIdentity( const IndexType size )
{
    _MatrixStorage::setDimension( size, size );

    SCAI_LOG_INFO( logger, "set identity values for " << size )

    HArrayUtils::setOrder( mIA, size, getContextPtr() );
    HArrayUtils::setOrder( mJA, size, getContextPtr() );

    HArrayUtils::setSameValue( mValues, size, ValueType( 1 ), getContextPtr() );
}

/* --------------------------------------------------------------------------- */
template<typename ValueType>
void COOStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    _MatrixStorage::setDimension( size, size );

    HArrayUtils::setOrder( mIA, size, getContextPtr() );
    HArrayUtils::setOrder( mJA, size, getContextPtr() );

    HArrayUtils::setArray( mValues, diagonal, common::BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setCOOData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    // check the sizes of the arrays

    SCAI_ASSERT_EQ_ERROR( numValues, ia.size(), "illegal size for ia" )
    SCAI_ASSERT_EQ_ERROR( numValues, ja.size() , "illegal size for ja" )
    SCAI_ASSERT_EQ_ERROR( numValues, values.size() , "illegal size for values" )

    _MatrixStorage::setDimension( numRows, numColumns );

    ContextPtr loc = getContextPtr();

    HArrayUtils::assign( mIA, ia, loc );
    HArrayUtils::assign( mJA, ja, loc );
    HArrayUtils::_assign( mValues, values, loc ); // supports type conversion

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode
#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "COOStorage.setCOOData" );
#endif

    // Note: no support for row indexes in COO format
    SCAI_LOG_INFO( logger, *this << ": set COO by arrays ia, ja, values" )

    verifySorting();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::purge()
{
    _MatrixStorage::setDimension( 0, 0 );

    mIA.purge();
    mJA.purge();
    mValues.purge();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate COO sparse matrix of size " << numRows << " x " << numColumns )

    _MatrixStorage::setDimension( numRows, numColumns );

    mIA.clear();
    mJA.clear();
    mValues.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "COOStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << getNumRows() << " x " << getNumColumns()
           << ", nnz = " << getNumValues() << " )" ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    IndexType pos = COOUtils::getValuePos( i, j, mIA, mJA, getContextPtr() );

    ValueType val = 0;

    if ( pos != invalidIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mValues.size(), "illegal value position for ( " << i << ", " << j << " )" );

        val = utilskernel::HArrayUtils::getVal<ValueType>( mValues, pos );
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    IndexType pos = COOUtils::getValuePos( i, j, mIA, mJA, getContextPtr() );

    if ( pos == invalidIndex )
    {
        COMMON_THROWEXCEPTION( "COO storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mValues.size(), "illegal value position for ( " << i << ", " << j << " )" );

    utilskernel::HArrayUtils::setVal( mValues, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mIA.prefetch( location );
    mJA.prefetch( location );
    mValues.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& COOStorage<ValueType>::getIA() const
{
    return mIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& COOStorage<ValueType>::getJA() const
{
    return mJA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& COOStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType COOStorage<ValueType>::getNumValues() const
{
    return mIA.size();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setDiagonal( const ValueType value )
{
    const IndexType M = getNumRows();
    const IndexType N = getNumColumns();

    SCAI_ASSERT_ERROR( COOUtils::hasDiagonal( M, N, mIA, mJA, getContextPtr() ),
                       "cannot set diagonal for COO, no all entries are available" )

    COOUtils::setDiagonal( mValues, value, M, N, mIA, mJA,  getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::UnaryOp::CONJ, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::scaleRows( const HArray<ValueType>& values )
{
    COOUtils::setRows( mValues, getNumRows(), getNumColumns(), mIA, mJA,
                       values, common::BinaryOp::MULT, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::scaleColumns( const HArray<ValueType>& values )
{
    COOUtils::setColumns( mValues, getNumRows(), getNumColumns(), mIA, mJA,
                          values, common::BinaryOp::MULT, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::wait() const
{
    mIA.wait();
    mJA.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::swap( COOStorage<ValueType>& other )
{
    // swap base class

    MatrixStorage<ValueType>::swap( other );

    // swap member varibles

    mIA.swap( other.mIA );
    mJA.swap( other.mJA );
    mValues.swap( other.mValues );

    // verify that row and column indexes in ia, ja are legal

#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "COOStorage.swap( ia, ja, values)" );
#endif

}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::splitUp( 
    IndexType& numRows, 
    IndexType& numColumns, 
    HArray<IndexType>& ia, 
    HArray<IndexType>& ja, 
    HArray<ValueType>& values )
{
    numRows    = getNumRows();
    numColumns = getNumColumns();

    ia = std::move( mIA );
    ja = std::move( mJA );
    values = std::move( mValues );

    // this storage should never be used afterwards, but make it consistent for safety

    _MatrixStorage::setDimension( 0, 0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> COOStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )
    // asum over the full array mValues
    return HArrayUtils::l1Norm( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> COOStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )
    ValueType res = HArrayUtils::dotProduct( mValues, mValues, this->getContextPtr() );
    return common::Math::sqrt( res );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> COOStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )
    return HArrayUtils::reduce( mValues, BinaryOp::ABS_MAX, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t COOStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( IndexType );
    memoryUsage += sizeof( IndexType ) * mIA.size();
    memoryUsage += sizeof( IndexType ) * mJA.size();
    memoryUsage += sizeof( ValueType ) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::incGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "incGEMV ( async = " << async << " ) , result += " << alpha << " * storage * x" )

    return COOUtils::gemv( result, 
                           getNumRows(), getNumColumns(), mIA, mJA, mValues, op, 
                           alpha, x, ValueType( 1 ), result, 
                           async, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const

{
    bool async = false; // synchronously execution, no SyncToken required
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    bool async = true;
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token, "NULL token not allowed for asynchronous execution gemv, alpha = " << alpha << ", beta = " << beta )
    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_REGION( "Storage.COO.gemv" )

    const IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    SCAI_LOG_INFO( logger,
                   "GEMV ( op = " << op << ", async = " << async
                   << " ), result [" << nTarget << "] = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::Constants::ZERO || getNumValues() == 0 )
    {
        // so we just have result = beta * y, will be done synchronously

        HArrayUtils::compute( result, beta, BinaryOp::MULT, y, getContextPtr() );

        if ( async )
        {
            return new tasking::NoSyncToken();
        }
        else
        {
            return NULL;
        }
    }

    // Due to COO format GEVM does not benefit of coupling all in one operation, so split it

    // Step 1: result = beta * y

    if ( beta == common::Constants::ZERO )
    {
        HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), getContextPtr() );
    }
    else
    {
        // Note: BinaryOp::MULT will deal with
        SCAI_ASSERT_EQUAL( y.size(), nTarget, "size mismatch y, beta = " << beta )
        HArrayUtils::compute( result, beta, BinaryOp::MULT, y, getContextPtr() );
    }

    return incGEMV( result, alpha, x, op, async );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const hmemo::HArray<OtherValueType>& values )
{
    IndexType numValues = ja.size();

    SCAI_LOG_DEBUG( logger, "set CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues )

    if ( ia.size() == numRows )
    {
        // offset array required
        hmemo::HArray<IndexType> offsets;
        IndexType total = CSRUtils::sizes2offsets( offsets, ia, getContextPtr() );
        SCAI_ASSERT_EQUAL( numValues, total, "sizes do not sum to number of values" );
        setCSRDataImpl( numRows, numColumns, offsets, ja, values );
        return;
    }

    _MatrixStorage::setDimension( numRows, numColumns );

    SCAI_REGION( "Storage.COO.buildCSR" )

    SCAI_ASSERT_EQUAL_DEBUG( numRows + 1, ia.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, values.size() )

    SCAI_LOG_DEBUG( logger,
                    "input csr data with " << numValues << "entries" )

    COOUtils::convertCSR2COO( mIA, ia, numValues, getContextPtr() );

    HArrayUtils::setArray( mJA, ja, common::BinaryOp::COPY, getContextPtr() );
    HArrayUtils::setArray( mValues, values, common::BinaryOp::COPY, getContextPtr() );

    verifySorting();
}

template<typename ValueType>
void COOStorage<ValueType>::setDiagonalV( const hmemo::HArray<ValueType>& diagonal )
{
    const IndexType M = getNumRows();
    const IndexType N = getNumColumns();

    SCAI_ASSERT_ERROR( COOUtils::hasDiagonal( M, N, mIA, mJA, getContextPtr() ), 
                       "cannot set diagonal for COO, not all diagonal entries available" )

    COOUtils::setDiagonalV( mValues, diagonal, M, N, mIA, mJA,  getContextPtr() );
}

template<typename ValueType>
void COOStorage<ValueType>::getDiagonal( hmemo::HArray<ValueType>& diagonal ) const
{
    COOUtils::getDiagonal( diagonal, getNumRows(), getNumColumns(), mIA, mJA, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.COO.getSparseRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    // resize the output arrays, invalidate old data before

    IndexType offset;    // first pos in COO data for row i
    IndexType n;         // number of entries in row i

    COOUtils::getRowPositions( offset, n, mIA, i, getContextPtr() );
 
    if ( n == 0 )
    {
        // no entries at all for row i

        jA.clear();
        values.clear();
        return;
    }
 
    // sparse data of row i is stored contiguously in cooJA and cooValues, so use setArraySection

    jA.resize( n );
    values.resize( n );

    HArrayUtils::setArraySection( jA, 0, 1, mJA, offset, 1, n, common::BinaryOp::COPY, getContextPtr() );
    HArrayUtils::setArraySection( values, 0, 1, mValues, offset, 1, n, common::BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getRow( hmemo::HArray<ValueType>& row, const IndexType i ) const
{
    SCAI_REGION( "Storage.COO.getRow" )

    HArray<IndexType> ja;
    HArray<ValueType> values;

    getSparseRow( ja, values, i );

    row.setSameValue( getNumColumns(), ValueType( 0 ) );

    HArrayUtils::scatter( row, ja, true, values, BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const BinaryOp op )
{
    SCAI_REGION( "Storage.COO.setRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    IndexType offset;
    IndexType n;
 
    COOUtils::getRowPositions( offset, n, mIA, i, getContextPtr() );

    HArray<IndexType> colIndexes;   // column indexes for row i

    colIndexes.resize( n );   // size must be set before as we deal with a section

    HArrayUtils::setArraySection( colIndexes, 0, 1, mJA, offset, 1, n, BinaryOp::COPY, getContextPtr() );

    HArray<ValueType> rowValues;    // contains the values of entries belonging to row i

    HArrayUtils::gather( rowValues, row, colIndexes, BinaryOp::COPY, getContextPtr() );

    HArrayUtils::setArraySection( mValues, offset, 1, rowValues, 0, 1, n, op, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& ia, hmemo::HArray<ValueType>& values, const IndexType j ) const
{   
    SCAI_REGION( "Storage.COO.getSparseCol" )
    
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )
    
    HArray<IndexType> pos;     // positions in the COO arrays with cooJA[pos] == j 

    COOUtils::getColumnPositions( pos, mJA, j, getContextPtr() );

    // column[ row ] = mValues[ pos ];

    HArrayUtils::gather( values, mValues, pos, BinaryOp::COPY, getContextPtr() );
    HArrayUtils::gather( ia, mIA, pos, BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.COO.getDenseCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, getNumRows(), colValues, rowIndexes, ValueType( 0 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const BinaryOp op )
{
    SCAI_REGION( "Storage.COO.setCol" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), getNumRows(), "column array to small for set" )

    HArray<IndexType> pos;    // positions in the ia, values array, are unique

    COOUtils::getColumnPositions( pos, mJA, j, getContextPtr() );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j

    HArrayUtils::gather( rowIndexes, mIA, pos, BinaryOp::COPY, getContextPtr() );

    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    HArrayUtils::gather( colValues, column, rowIndexes, BinaryOp::COPY, getContextPtr() );

    HArrayUtils::scatter( mValues, pos, true, colValues, op, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::scale( const ValueType value )
{
    // multiply value with each entry of mValues
    utilskernel::HArrayUtils::setScalar( mValues, value, BinaryOp::MULT, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.COO.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), getNumColumns() )              // matrix must be square

    COOUtils::jacobi( solution, omega, oldSolution, rhs,
                      mIA, mJA, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.COO.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    SCAI_ASSERT_EQ_ERROR( getNumRows(), localSolution.size(), "array localSolution has illegal size" )
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), oldHaloSolution.size(), "array old halo solution has illegal size" )

    COOUtils::jacobiHalo( localSolution, omega, localDiagonal, oldHaloSolution,
                          mIA, mJA, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::buildCSR(
    hmemo::HArray<IndexType>& csrIA,
    hmemo::HArray<IndexType>* csrJA,
    hmemo::HArray<OtherValueType>* csrValues,
    const hmemo::ContextPtr preferredLoc ) const
{
    SCAI_REGION( "Storage.COO.buildCSR" )

    if ( csrJA == NULL || csrValues == NULL )
    {
        // number of entries per row, count with buckets for each row

        SCAI_LOG_INFO( logger, "build CSR sizes from " << *this )

        utilskernel::HArrayUtils::bucketCount( csrIA, mIA, getNumRows(), preferredLoc );
        return;
    }

    SCAI_LOG_INFO( logger, "build CSR data from " << *this )

    COOUtils::convertCOO2CSR( csrIA, mIA, getNumRows(), preferredLoc );

    SCAI_ASSERT_EQ_DEBUG( getNumRows() + 1, csrIA.size(), "serious mismatch, should not happen" )

    // CSR array ja, values can be directly copied

    utilskernel::HArrayUtils::setArray( *csrJA, mJA, BinaryOp::COPY, preferredLoc );
    utilskernel::HArrayUtils::setArray( *csrValues, mValues, BinaryOp::COPY, preferredLoc );

    // Note: sort is stable, so diagonal values remain first in each row
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const ValueType beta,
    const MatrixStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", with A = " << a << ", B = " << b )

    SCAI_REGION( "Storage.COO.plusMatrix" )

    if ( a.getFormat() != Format::COO )
    {
        matrixPlusMatrix( alpha, convert<COOStorage<ValueType>>( a ), beta, b );
    } 
    else if ( b.getFormat() != Format::COO )
    {
        matrixPlusMatrix( alpha, a, beta, convert<COOStorage<ValueType>>( b ) );
    }
    else
    {
        // a and b have the right format so we just cast it correctly

        matrixPlusMatrixImpl( alpha, static_cast<const COOStorage<ValueType>&>( a ), 
                              beta, static_cast<const COOStorage<ValueType>&>( b ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::matrixPlusMatrixImpl(
    const ValueType alpha,
    const COOStorage<ValueType>& a,
    const ValueType beta,
    const COOStorage<ValueType>& b )
{
    if ( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary
        COOStorage<ValueType> tmp;
        tmp.matrixPlusMatrixImpl( alpha, a, beta, b );
        SCAI_LOG_DEBUG( logger, "swap this = " << *this << " and tmp = " << tmp )
        swap( tmp ); // safe as tmp will be destroyed afterwards
        return;
    }

    SCAI_REGION( "Storage.COO.plusImpl" )

    // for COO we do not require that a and b have same size

    _MatrixStorage::setDimension( common::Math::max( a.getNumRows(), b.getNumRows() ),
                                  common::Math::max( a.getNumColumns(), b.getNumColumns() ) );

    const IndexType nnz1 = a.getNumValues();
    const IndexType nnz2 = b.getNumValues();

    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " * A + " << beta << " * B, with " << "A = " << a << ", B = " << b << ", all are COO" )

    // concatenae the COO data of both storages

    mIA.resize( nnz1 + nnz2 );
    mJA.resize( nnz1 + nnz2 );
    mValues.resize( nnz1 + nnz2 );

    HArrayUtils::setArraySection( mIA, 0, 1, a.getIA(), 0, 1, nnz1 );
    HArrayUtils::setArraySection( mIA, nnz1, 1, b.getIA(), 0, 1, nnz2 );
    HArrayUtils::setArraySection( mJA, 0, 1, a.getJA(), 0, 1, nnz1 );
    HArrayUtils::setArraySection( mJA, nnz1, 1, b.getJA(), 0, 1, nnz2 );

    HArray<ValueType> tmp;

    HArrayUtils::binaryOpScalar( tmp, a.getValues(), alpha, BinaryOp::MULT, false );
    HArrayUtils::setArraySection( mValues, 0, 1, tmp, 0, 1, nnz1 );
    HArrayUtils::binaryOpScalar( tmp, b.getValues(), beta, BinaryOp::MULT, false );
    HArrayUtils::setArraySection( mValues, nnz1, 1, tmp, 0, 1, nnz2 );

    // sort it and make it unique

    COOUtils::normalize( mIA, mJA, mValues, BinaryOp::ADD, getContextPtr() );

    SCAI_LOG_INFO( logger, "COO matrix add: nnz = " << getNumValues() << " from " << nnz1 << " + " << nnz2 )
}

/* ========================================================================= */
/*  Halo stuff                                                               */
/* ========================================================================= */

template<typename ValueType>
void COOStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{
    halo.halo2Global( mJA );
    _MatrixStorage::setDimension( getNumRows(), globalNumColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::copy() const
{
    return new COOStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    std::unique_ptr<COOStorage<ValueType> > storage( new COOStorage<ValueType>( getContextPtr() ) );
    storage->allocate( numRows, numColumns );
    return storage.release();
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

template<typename ValueType>
std::string COOStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "COOStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* COOStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
const char* COOStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType COOStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::COO, common::getScalarType<ValueType>() );
}

template<typename ValueType>
MatrixStorageCreateKeyType COOStorage<ValueType>::getCreateValue() const
{
    return createValue();
}

template<typename ValueType>
_MatrixStorage* COOStorage<ValueType>::create()
{
    return new COOStorage<ValueType>();
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */


SCAI_COMMON_INST_CLASS( COOStorage, SCAI_NUMERIC_TYPES_HOST )

#define COO_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                   \
    template void COOStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,     \
            hmemo::HArray<OtherValueType>* values,const hmemo::ContextPtr ) const;                           \
    template void COOStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other );           \
    template void COOStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType,                   \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                \
            const hmemo::HArray<OtherValueType>& );                          

#define COO_STORAGE_INST_LVL1( ValueType )                                                                   \
    SCAI_COMMON_LOOP_LVL2( ValueType, COO_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( COO_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef COO_STORAGE_INST_LVL2
#undef COO_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
