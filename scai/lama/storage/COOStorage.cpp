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
 * @author Lauretta Schubert
 * @date 25.05.2011
 */

// hpp
#include <scai/lama/storage/COOStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/loop.hpp>

// sqrt for all value types
#include <cmath>

using namespace scai::hmemo;

namespace scai
{

using common::unique_ptr;
using common::shared_ptr;

using tasking::SyncToken;

using utilskernel::UtilKernelTrait;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;

using sparsekernel::COOKernelTrait;
using sparsekernel::CSRKernelTrait;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, COOStorage<ValueType>::logger, "MatrixStorage.COOStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( const IndexType numRows, const IndexType numColumns )
    :

    CRTPMatrixStorage<COOStorage<ValueType>, ValueType>( numRows, numColumns ), mNumValues( 0 )
{
    SCAI_LOG_DEBUG( logger, "COOStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values )

    : CRTPMatrixStorage<COOStorage<ValueType>, ValueType>()
{
    // all array must have the same size
    IndexType numValues = ia.size();
    setCOOData( numRows, numColumns, numValues, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( const COOStorage<ValueType>& other )

    : CRTPMatrixStorage<COOStorage<ValueType>, ValueType>( 0, 0 )
{
    // ToDo: copy of same storage format should be more efficient
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage()
    : CRTPMatrixStorage<COOStorage<ValueType>, ValueType>( 0, 0 ), mNumValues( 0 )
{
    SCAI_LOG_DEBUG( logger, "COOStorage, matrix is 0 x 0." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format::MatrixStorageFormat COOStorage<ValueType>::getFormat() const
{
    return Format::COO;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "COOStorage " << mNumRows << " x " << mNumColumns << ", #values = " << mNumValues << endl;
    ContextPtr host = Context::getHostPtr();
    ReadAccess<IndexType> ia( mIA, host );
    ReadAccess<IndexType> ja( mJA, host );
    ReadAccess<ValueType> values( mValues, host );

    for ( IndexType i = 0; i < mNumValues; i++ )
    {
        stream << "@[ " << ia[i] << ", " << ja[i] << " ] = " << values[i] << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool COOStorage<ValueType>::checkDiagonalProperty() const
{
    bool diagonalProperty = true;

    if ( mNumRows != mNumColumns )
    {
        diagonalProperty = false;
    }
    else if ( mNumRows == 0 )
    {
        // zero sized matrix
        diagonalProperty = true;
    }
    else if ( mIA.size() == 0 )
    {
        diagonalProperty = false;
    }
    else
    {
        diagonalProperty = true; // intialization for reduction
        static LAMAKernel<COOKernelTrait::hasDiagonalProperty> hasDiagonalProperty;
        ContextPtr loc = this->getContextPtr();
        hasDiagonalProperty.getSupportedContext( loc );
        ReadAccess<IndexType> ia( mIA, loc );
        ReadAccess<IndexType> ja( mJA, loc );
        diagonalProperty = hasDiagonalProperty[loc]( ia.get(), ja.get(), mNumRows );
    }

    SCAI_LOG_INFO( logger, *this << ": checkDiagonalProperty -> " << diagonalProperty )
    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;
    mNumValues = 0;
    mIA.clear();
    mJA.clear();
    mValues.clear();
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::check( const char* msg ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mIA.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mJA.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )
    // check row indexes in IA and column indexes in JA
    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
        ContextPtr loc = getContextPtr();
        validIndexes.getSupportedContext( loc );  // find location where routine is available
        ReadAccess<IndexType> rJA( mJA, loc );
        ReadAccess<IndexType> rIA( mIA, loc );
        SCAI_CONTEXT_ACCESS( loc )
        bool okayIA = validIndexes[ loc ]( rIA.get(), mNumValues, mNumRows );
        SCAI_ASSERT_ERROR( okayIA,  *this << " @ " << msg << ": illegel indexes in IA" )
        bool okayJA = validIndexes[ loc ]( rJA.get(), mNumValues, mNumColumns );
        SCAI_ASSERT_ERROR( okayJA, *this << " @ " << msg << ": illegel indexes in JA" )
        SCAI_LOG_INFO( logger, "check, msg = " << msg << ", okayIA = " << okayIA << ", okayJA = " << okayJA )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_INFO( logger, "set identity values for " << size )
    mNumRows = size;
    mNumColumns = size;
    mNumValues = mNumRows;
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    ContextPtr loc = this->getContextPtr();
    setOrder.getSupportedContext( loc, setVal );  // setOrder, setVal must be available at loc
    WriteOnlyAccess<IndexType> ia( mIA, loc, mNumValues );
    WriteOnlyAccess<IndexType> ja( mJA, loc, mNumValues );
    WriteOnlyAccess<ValueType> values( mValues, loc, mNumValues );
    SCAI_CONTEXT_ACCESS( loc )
    setOrder[loc]( ia.get(), mNumValues );
    setOrder[loc]( ja.get(), mNumValues );
    setVal[loc]( values.get(), mNumValues, ValueType( 1 ), common::binary::COPY );
    mDiagonalProperty = true;
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
    SCAI_ASSERT_EQUAL_ERROR( numValues, ia.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, values.size() )
    _MatrixStorage::setDimension( numRows, numColumns );
    mNumValues = numValues;
    ContextPtr loc = getContextPtr();
    HArrayUtils::assign( mIA, ia, loc );
    HArrayUtils::assign( mJA, ja, loc );
    HArrayUtils::assign( mValues, values, loc ); // supports type conversion
    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode
#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "COOStorage.setCOOData" );
#endif
    mDiagonalProperty = checkDiagonalProperty();
    // Note: no support for row indexes in COO format
    SCAI_LOG_INFO( logger, *this << ": set COO by arrays ia, ja, values" )
}

/* --------------------------------------------------------------------------- */



/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::~COOStorage()
{
    SCAI_LOG_DEBUG( logger, "~COOStorage for matrix " << mNumRows << " x " << mNumColumns )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::purge()
{
    mNumColumns = 0;
    mNumRows = 0;
    mNumValues = 0;
    mIA.purge();
    mJA.purge();
    mValues.purge();
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate COO sparse matrix of size " << numRows << " x " << numColumns )
    clear(); // all variables are set for a zero-sized matrix
    mNumRows = numRows;
    mNumColumns = numColumns;
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "COOStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << mNumRows << " x " << mNumColumns
           << ", nnz = " << mNumValues << " )" ;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<COOKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get(), mNumValues );

    ValueType val = 0;

    if ( pos != nIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumValues, "illegal value position for ( " << i << ", " << j << " )" );

        val = utilskernel::HArrayUtils::getVal<ValueType>( mValues, pos );
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<COOKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get(), mNumValues );

    if ( pos == nIndex )
    {
        COMMON_THROWEXCEPTION( "COO storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumValues, "illegal value position for ( " << i << ", " << j << " )" );

    utilskernel::HArrayUtils::setValImpl( mValues, pos, val, op );
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
    return mNumValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
    ContextPtr loc = Context::getHostPtr();
    WriteAccess<ValueType> wValues( mValues, loc );
    ReadAccess<IndexType> rJa( mJA, loc );
    ReadAccess<IndexType> rIa( mIA, loc );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wValues[i] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::unary::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::scaleImpl( const HArray<OtherType>& values )
{
    static LAMAKernel<COOKernelTrait::scaleRows<ValueType, OtherType> > scaleRows;
    ContextPtr loc = this->getContextPtr();
    scaleRows.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<OtherType> rValues( values, loc );
    WriteAccess<ValueType> wValues( mValues, loc );  // update
    ReadAccess<IndexType> rIa( mIA, loc );
    scaleRows[loc]( wValues.get(), rValues.get(), rIa.get(), mNumValues );
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
void COOStorage<ValueType>::swap( _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( getFormat(), other.getFormat(), "swap only for same storage format" )
    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(), "swap only for same value type" )

    // only in debug mode use the more expensive dynamic cast for verification

    SCAI_ASSERT_DEBUG( dynamic_cast<COOStorage<ValueType>* >( &other ), "illegal storage to swap" )

    swapImpl( reinterpret_cast<COOStorage<ValueType>& >( other ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::swapImpl( COOStorage<ValueType>& other )
{
    std::swap( mNumValues, other.mNumValues );
    mIA.swap( other.mIA );
    mJA.swap( other.mJA );
    mValues.swap( other.mValues );

    MatrixStorage<ValueType>::swapMS( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::swap( HArray<IndexType>& ia, HArray<IndexType>& ja, HArray<ValueType>& values )
{
    IndexType numValues = ia.size();

    SCAI_ASSERT_EQUAL( numValues, ja.size(), "mismatch of coo IA and JA array" );
    SCAI_ASSERT_EQUAL( numValues, values.size(), "mismatch of coo IA and Values array" );

    mIA.swap( ia );
    mJA.swap( ja );
    mValues.swap( values );

    // guarantee consistency of the new array

    mNumValues = numValues;
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )
    // asum over the full array mValues
    return HArrayUtils::asum( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )
    ValueType res = HArrayUtils::dotProduct( mValues, mValues, this->getContextPtr() );
    return common::Math::sqrt( res );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )
    return HArrayUtils::reduce( mValues, common::binary::ABS_MAX, this->getContextPtr() );
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
    bool async ) const
{
    SCAI_LOG_INFO( logger, "incGEMV ( async = " << async << " ) , result += " << alpha << " * storage * x" )
    static LAMAKernel<COOKernelTrait::normalGEMV<ValueType> > normalGEMV;
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    common::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );
    // use general kernel, might change
    ValueType beta = 1;
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, wResult.get(),
                     mNumRows, mNumValues,
                     cooIA.get(), cooJA.get(), cooValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( cooValues.releaseDelayed() );
        syncToken->pushRoutine( cooIA.releaseDelayed() );
        syncToken->pushRoutine( cooJA.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::matrixTimesVector(

    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const

{
    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.COO.matrixTimesVector" )

    if ( alpha == common::constants::ZERO || mNumValues == 0 )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::binaryOpScalar1( result, beta, y, common::binary::MULT, this->getContextPtr() );
        return;
    }

    ContextPtr loc = this->getContextPtr();

    // Due to COO format GEMV does not benefit of coupling all in one operation, so split it

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumRows );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
        HArrayUtils::binaryOpScalar1( result, beta, y, common::binary::MULT, this->getContextPtr() );
    }

    bool async = false;
    SyncToken* token = incGEMV( result, alpha, x, async );
    SCAI_ASSERT( token == NULL, "syncrhonous execution cannot have token" )
}


/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::incGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "incGEVM ( async = " << async << " ) , result += " << alpha << " * x * storage" )
    static LAMAKernel<COOKernelTrait::normalGEVM<ValueType> > normalGEVM;
    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );
    common::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc, mNumColumns );
    // use general kernel, might change
    ValueType beta = 1;
    normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, wResult.get(),
                     mNumColumns, mNumValues,
                     cooIA.get(), cooJA.get(), cooValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( cooValues.releaseDelayed() );
        syncToken->pushRoutine( cooIA.releaseDelayed() );
        syncToken->pushRoutine( cooJA.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.COO.VectorTimesMatrix" )
    ContextPtr loc = this->getContextPtr();

    // Due to COO format GEVM does not benefit of coupling all in one operation, so split it

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumColumns );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
        HArrayUtils::binaryOpScalar1( result, beta, y, common::binary::MULT, loc );
    }

    // Step 2: result = alpha * x * this + 1 * result
    bool async = false;
    SyncToken* token = incGEVM( result, alpha, x, async );
    SCAI_ASSERT( NULL == token, "syncrhonous execution has no token" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVectorAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.COO.matrixTimesVectorAsync" )
    ContextPtr loc = this->getContextPtr();

    if ( alpha == common::constants::ZERO || mNumValues == 0 )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::binaryOpScalar1( result, beta, y, common::binary::MULT, loc );
        return new tasking::NoSyncToken();
    }

    // Due to COO format GEVM does not benefit of coupling all in one operation, so split it

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumRows );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        // Note: binary::MULT will deal with
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
        HArrayUtils::binaryOpScalar1( result, beta, y, common::binary::MULT, loc );
    }

    bool async = true;
    SyncToken* token = incGEMV( result, alpha, x, async );
    SCAI_ASSERT( token, "asyncrhonous execution cannot have NULL token" )
    return token;
}

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const hmemo::HArray<OtherValueType>& values,
    const hmemo::ContextPtr prefLoc )
{
    SCAI_LOG_DEBUG( logger, "set CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues )

    if ( ia.size() == numRows )
    {
        // offset array required
        hmemo::HArray<IndexType> offsets;
        IndexType total = _MatrixStorage::sizes2offsets( offsets, ia, prefLoc );
        SCAI_ASSERT_EQUAL( numValues, total, "sizes do not sum to number of values" );
        setCSRDataImpl( numRows, numColumns, numValues, offsets, ja, values, prefLoc );
        return;
    }

    SCAI_REGION( "Storage.COO.buildCSR" )

    SCAI_ASSERT_EQUAL_DEBUG( numRows + 1, ia.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, values.size() )
    hmemo::ContextPtr loc = prefLoc;
    // ReadAccess<IndexType> csrJA( ja, loc );
    // ReadAccess<OtherValueType> csrValues( values, loc );
    mNumRows = numRows;
    mNumColumns = numColumns;
    // check if input csr data has the diagonal property and inherit it
    int numDiagonals = std::min( numRows, numColumns );
    {
        SCAI_LOG_DEBUG( logger,
                        "check CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues << " for diagonal property, #diagonals = " << numDiagonals )
        static utilskernel::LAMAKernel<sparsekernel::CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;
        hmemo::ContextPtr loc = this->getContextPtr();
        hasDiagonalProperty.getSupportedContext( loc );
        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::ReadAccess<IndexType> csrJA( ja, loc );
        SCAI_CONTEXT_ACCESS( loc )
        mDiagonalProperty = hasDiagonalProperty[loc] ( numDiagonals, csrIA.get(), csrJA.get() );
    }

    if ( !mDiagonalProperty )
    {
        numDiagonals = 0; // do not store diagonal data at the beginning in COO data
    }

    mNumValues = numValues;
    SCAI_LOG_DEBUG( logger,
                    "input csr data with " << mNumValues << "entries,  has diagonal property = " << mDiagonalProperty )
    {
        static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::offsets2ia> offsets2ia;
        hmemo::ContextPtr loc = this->getContextPtr();
        offsets2ia.getSupportedContext( loc );
        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::WriteOnlyAccess<IndexType> cooIA( mIA, loc, mNumValues );
        SCAI_CONTEXT_ACCESS( loc )
        offsets2ia[loc]( cooIA.get(), mNumValues, csrIA.get(), mNumRows, numDiagonals );
    }
    {
        static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::setCSRData<IndexType, IndexType> > setCSRData;
        hmemo::ContextPtr loc = this->getContextPtr();   // preferred location
        setCSRData.getSupportedContext( loc );    // supported location
        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::ReadAccess<IndexType> csrJA( ja, loc );
        hmemo::WriteOnlyAccess<IndexType> cooJA( mJA, loc, mNumValues );
        SCAI_CONTEXT_ACCESS( loc )
        setCSRData[loc]( cooJA.get(), csrJA.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }
    {
        static utilskernel::LAMAKernel<sparsekernel::COOKernelTrait::setCSRData<ValueType, OtherValueType> > setCSRData;
        hmemo::ContextPtr loc = this->getContextPtr();   // preferred location
        setCSRData.getSupportedContext( loc );    // supported location
        hmemo::ReadAccess<IndexType> csrIA( ia, loc );
        hmemo::ReadAccess<OtherValueType> csrValues( values, loc );
        hmemo::WriteOnlyAccess<ValueType> cooValues( mValues, loc, mNumValues );
        SCAI_CONTEXT_ACCESS( loc )
        setCSRData[loc]( cooValues.get(), csrValues.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }
}

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::setDIADataImpl(
    const IndexType /*numRows*/,
    const IndexType /*numColumns*/,
    const IndexType /*numDiagonals*/,
    const HArray<IndexType>& /*offsets*/,
    const HArray<OtherValueType>& /*values*/,
    const ContextPtr /*prefLoc*/ )
{
    COMMON_THROWEXCEPTION( "not yet implemeted" )
}

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherType>& diagonal )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<ValueType, OtherType> > set;
    hmemo::ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    hmemo::ReadAccess<OtherType> rDiagonal( diagonal, loc );
    hmemo::WriteAccess<ValueType> wValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    // diagonal elements are the first entries of mValues
    set[loc]( wValues.get(), rDiagonal.get(), numDiagonalElements, common::binary::COPY );
}

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherType>& diagonal ) const
{
    // diagional[0:numDiagonalElements] = mValues[0:numDiagonalElements]
    // Note: using HArrayUtils::setArray not possible, as we only need part of mValues
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
    static utilskernel::LAMAKernel<utilskernel::UtilKernelTrait::set<OtherType, ValueType> > set;
    hmemo::ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    hmemo::WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    hmemo::ReadAccess<ValueType> rValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    // diagonal elements are the first entries of mValues
    set[loc]( wDiagonal.get(), rValues.get(), numDiagonalElements, common::binary::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::_HArray& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.COO.getSparseRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )

    // resize the output arrays, invalidate old data before

    static LAMAKernel<COOKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = this->getContextPtr();
    getValuePosRow.getSupportedContext( loc );

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wColIndexes( jA, loc, mNumColumns );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumColumns );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosRow[loc]( wColIndexes.get(), wValuePos.get(), i,
                                             rIA.get(), mNumColumns, rJA.get(), mNumValues );

        wColIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    values.clear();
    values.resize( valuePos.size());

    HArrayUtils::gather( values, mValues, valuePos, common::binary::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getRowImpl( hmemo::HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_REGION( "Storage.COO.getRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )

    static LAMAKernel<COOKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = this->getContextPtr();

    getValuePosRow.getSupportedContext( loc );

    HArray<IndexType> colIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> rowValues;    // contains the values of entries belonging to row i

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wColIndexes( colIndexes, loc, mNumColumns );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumColumns );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosRow[loc]( wColIndexes.get(), wValuePos.get(), i,
                                             rIA.get(), mNumColumns, rJA.get(), mNumValues );

        wColIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    row.init( ValueType( 0 ), mNumColumns );

    // row[ colIndexes ] = mValues[ pos ];

    HArrayUtils::gatherImpl( rowValues, mValues, valuePos, common::binary::COPY, loc );
    HArrayUtils::scatterImpl( row, colIndexes, true, rowValues, common::binary::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::setRowImpl( const HArray<OtherType>& row, const IndexType i,
                                        const common::binary::BinaryOp op )
{
    SCAI_REGION( "Storage.COO.setRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), mNumColumns, "row array to small for set" )

    // get sparse vector with column indexes and positions

    static LAMAKernel<COOKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = this->getContextPtr();

    getValuePosRow.getSupportedContext( loc );

    HArray<IndexType> colIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> rowValues;    // contains the values of entries belonging to row i

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wColIndexes( colIndexes, loc, mNumColumns );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumColumns );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosRow[loc]( wColIndexes.get(), wValuePos.get(), i,
                                             rIA.get(), mNumColumns, rJA.get(), mNumValues );

        wColIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    // mValues[pos] = row[ colIndexes ]

    HArrayUtils::gatherImpl( rowValues, row, colIndexes, common::binary::COPY, loc );
    HArrayUtils::scatterImpl( mValues, valuePos, true, rowValues, op, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::_HArray& values, const IndexType j ) const
{   
    SCAI_REGION( "Storage.COO.getSparseCol" )
    
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "col index out of range" )
    
    static LAMAKernel<COOKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( iA, loc, mNumRows );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumRows );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), mNumRows, rJA.get(), mNumValues );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    // column[ row ] = mValues[ pos ];

    HArrayUtils::gather( values, mValues, valuePos, common::binary::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::getColumn( _HArray& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.COO.getDenseCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, mNumRows, colValues, rowIndexes );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::setColumnImpl( const HArray<OtherType>& column, const IndexType j,
        const common::binary::BinaryOp op )
{
    SCAI_REGION( "Storage.COO.setCol" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), mNumRows, "column array to small for set" )

    static LAMAKernel<COOKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array, are unique
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, mNumRows );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumRows );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), mNumRows, rJA.get(), mNumValues );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    //  mValues[ pos ] op= column[ rowIndexes ]

    HArrayUtils::gatherImpl( colValues, column, rowIndexes, common::binary::COPY, loc );
    HArrayUtils::scatterImpl( mValues, valuePos, true, colValues, op, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::scaleImpl( const ValueType value )
{
    // multiply value with each entry of mValues
    utilskernel::HArrayUtils::setScalar( mValues, value, common::binary::MULT, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::vectorTimesMatrixAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.COO.vectorTimesMatrixAsync" )
    ContextPtr loc = this->getContextPtr();

    // Due to COO format GEVM does not benefit of coupling all in one operation, so split it

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumColumns );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
        HArrayUtils::binaryOpScalar1( result, beta, y, common::binary::MULT, loc );
    }

    // Step 2: result = alpha * x * this + 1 * result
    bool async = true;
    SyncToken* token = incGEVM( result, alpha, x, async );
    SCAI_ASSERT( token, "asyncrhonous execution cannot have NULL token" )
    return token;
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
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )              // matrix must be square
    static LAMAKernel<COOKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );
    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );
    // Due to diagonal property there is no advantage by taking row indexes
    SCAI_CONTEXT_ACCESS( loc )
    jacobi[loc]( wSolution.get(), mNumValues, cooIA.get(), cooJA.get(), cooValues.get(),
                 rOldSolution.get(), rRhs.get(), omega, mNumRows );
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

        utilskernel::HArrayUtils::bucketCount( csrIA, mIA, mNumRows, preferredLoc );
        return;
    }

    SCAI_LOG_INFO( logger, "build CSR data from " << *this )

    HArray<IndexType> perm;  // help array for resorting the values

    utilskernel::HArrayUtils::bucketSort( csrIA, perm, mIA, mNumRows );

    SCAI_ASSERT_EQ_DEBUG( mNumRows + 1, csrIA.size(), "serious mismatch, should not happen" )
    SCAI_ASSERT_EQ_ERROR( perm.size(), mIA.size(), "Illegal entries in mIA of COO storage" )

    // CSR array ja, values are the COO arrays resorted

    utilskernel::HArrayUtils::gatherImpl( *csrJA, mJA, perm, common::binary::COPY, preferredLoc );
    utilskernel::HArrayUtils::gatherImpl( *csrValues, mValues, perm, common::binary::COPY, preferredLoc );

    // Note: sort is stable, so diagonal values remain first in each row
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::copy() const
{
    return new COOStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::newMatrixStorage() const
{
    common::unique_ptr<COOStorage<ValueType> > storage( new COOStorage<ValueType>() );
    storage->setContextPtr( this->getContextPtr() );
    return storage.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* COOStorage<ValueType>::create()
{
    return new COOStorage<ValueType>();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType COOStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::COO, common::getScalarType<ValueType>() );
}

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

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( COOStorage, SCAI_NUMERIC_TYPES_HOST )

#define COO_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                 \
    template void COOStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
            hmemo::HArray<OtherValueType>* values,const hmemo::ContextPtr ) const;                                         \
    template void COOStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                              \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                               \
    template void COOStorage<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;              \
    template void COOStorage<ValueType>::setRowImpl( const hmemo::HArray<OtherValueType>&, const IndexType,                \
            const common::binary::BinaryOp );                          \
    template void COOStorage<ValueType>::setColumnImpl( const hmemo::HArray<OtherValueType>&, const IndexType,             \
            const common::binary::BinaryOp );                       \
    template void COOStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                          \
    template void COOStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                          \
    template void COOStorage<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                                \
    template void COOStorage<ValueType>::setDIADataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );

#define COO_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, COO_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( COO_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef COO_STORAGE_INST_LVL2
#undef COO_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
