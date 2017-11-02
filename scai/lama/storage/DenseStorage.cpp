/**
 * @file DenseStorage.cpp
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
 * @brief Instantiation for template class DenseStorage.
 * @author Thomas Brandes, Michael Drost
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/DIAStorage.hpp>

// local libraries
#include <scai/lama/mepr/DenseStorageWrapper.hpp>

// internal scai libraries
#include <scai/sparsekernel/DenseKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>
#include <cmath>


using std::abs;
// so we can use abs for float and double and own abs for Complex

using std::shared_ptr;

namespace scai
{

using common::TypeTraits;

using namespace hmemo;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using utilskernel::UtilKernelTrait;

using sparsekernel::DenseKernelTrait;
using sparsekernel::CSRKernelTrait;


namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseStorage<ValueType>::logger,
                              "MatrixStorage.DenseStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
HArray<ValueType>& DenseStorage<ValueType>::getData()
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& DenseStorage<ValueType>::getData() const
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DenseStorage<ValueType>::getNumValues() const
{
    static LAMAKernel<DenseKernelTrait::nonZeroValues<ValueType> > nonZeroValues;
    ContextPtr loc = this->getContextPtr();
    nonZeroValues.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> values( mData, loc );
    IndexType count = nonZeroValues[loc]( values.get(), mNumRows, mNumColumns, MatrixStorage<ValueType>::mEpsilon );
    SCAI_LOG_INFO( logger, *this << ": #non-zero values = " << count )
    return count;
}

/* --------------------------------------------------------------------------- */

#ifdef SCAI_ASSERT_LEVEL_OFF
template<typename ValueType>
void DenseStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void DenseStorage<ValueType>::check( const char* /* msg */ ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows * mNumColumns, mData.size() )
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    static LAMAKernel<DenseKernelTrait::setDiagonalValue<ValueType> > setDiagonalValue;
    ContextPtr loc = this->getContextPtr();
    setDiagonalValue.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wData( mData, loc ); // use existing data
    setDiagonalValue[loc]( wData.get(), mNumRows, mNumColumns, value );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::_HArray& values, const IndexType i ) const
{
    // ToDo: avoid temporary array row by new version buildSparseArray with offs and inc argument

    HArray<ValueType> row;
    getRowImpl( row, i );
    HArrayUtils::buildSparseArray( values, jA, row, mContext );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::_HArray& values, const IndexType j ) const
{
    HArray<ValueType> col;
    getColumn( col, j );
    HArrayUtils::buildSparseArray( values, iA, col, mContext );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorage<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType rowIndex ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rowIndex, mNumRows, "row index out of range" )

    row.clear();                    // make all data invalid
    row.resize( mNumColumns );      // resize it

    // inc = denseindex( i, j + 1, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( rowIndex, 0, numRows, numColumns )

    const IndexType inc   = 1;
    const IndexType first = rowIndex * mNumColumns;

    HArrayUtils::setArraySection( row, 0, 1,             // row (:)
                                  mData, first, inc,
                                  mNumColumns, common::binary::COPY, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorage<ValueType>::setRowImpl( const HArray<OtherType>& row, const IndexType rowIndex,
        const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( rowIndex, mNumRows, "row index out of range" )

    // inc = denseindex( i, j + 1, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( i, 0, numRows, numColumns )

    const IndexType inc   = 1;
    const IndexType first = rowIndex * mNumColumns;

    HArrayUtils::setArraySection( mData, first, inc,
                                  row, 0, 1,
                                  mNumColumns,
                                  op,
                                  this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getColumn( _HArray& column, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    column.clear();                 // make all data invalid
    column.resize( mNumRows );      // resize it

    // inc = denseindex( i + 1, j, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( 0, j, numRows, numColumns )

    const IndexType inc   = mNumColumns;
    const IndexType first = j;

    HArrayUtils::setArraySection( column, 0, 1,
                                  mData, first, inc,
                                  mNumRows,
                                  common::binary::COPY,
                                  this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorage<ValueType>::setColumnImpl( const HArray<OtherType>& column, const IndexType colIndex,
        const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( colIndex, mNumColumns, "column index out of range" )

    // inc = denseindex( i + 1, j, numRows, numColumns ) - denseindex( i, j, numRows, numColumns )
    // first = denseindex( 0, j, numRows, numColumns )

    const IndexType inc   = mNumColumns;
    const IndexType first = colIndex;

    HArrayUtils::setArraySection( mData, first, inc ,
                                  column, 0, 1,
                                  mNumRows,
                                  op,
                                  this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorage<ValueType>::getDiagonalImpl( HArray<OtherType>& diagonal ) const
{
    IndexType numDiagonalValues = common::Math::min( mNumColumns, mNumRows );

    diagonal.clear();                       // make all data invalid
    diagonal.resize( numDiagonalValues );   // resize it

    // inc = denseindex( i + 1, i + 1, numRows, numColumns ) - denseindex( i, i, numRows, numColumns )
    // first = denseindex( 0, 0, numRows, numColumns )

    const IndexType inc   = mNumRows + 1;
    const IndexType first = 0;

    HArrayUtils::setArraySection( diagonal, 0, 1,
                                  mData, first, inc,
                                  numDiagonalValues,
                                  common::binary::COPY,
                                  mData.getValidContext() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorage<ValueType>::setDiagonalImpl( const HArray<OtherType>& diagonal )
{
    IndexType numDiagonalValues = common::Math::min( mNumColumns, mNumRows );

    SCAI_ASSERT_GE_DEBUG( diagonal.size(), numDiagonalValues, "diagonal array has insufficient size" )

    // inc = denseindex( i + 1, i + 1, numRows, numColumns ) - denseindex( i, i, numRows, numColumns )
    // first = denseindex( 0, 0, numRows, numColumns )

    const IndexType inc   = mNumRows + 1;
    const IndexType first = 0;

    HArrayUtils::setArraySection( mData, first, inc,
                                  diagonal, 0, 1,
                                  numDiagonalValues,
                                  common::binary::COPY,
                                  mData.getValidContext() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::scaleImpl( const ValueType value )
{
    // not used here: HArrayUtils::scale( mData, value, this->getContextPtr() )
    // reasoning:     workload distribution would not fit to distribution of rows
    static LAMAKernel<DenseKernelTrait::setValue<ValueType> > setValue;
    ContextPtr loc = this->getContextPtr();
    setValue.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wData( mData, loc );
    setValue[loc]( wData.get(), mNumRows, mNumColumns, value, common::binary::MULT );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mData, mData, common::unary::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorage<ValueType>::scaleImpl( const HArray<OtherType>& values )
{
    static LAMAKernel<DenseKernelTrait::scaleRows<ValueType, OtherType> > scaleRows;
    ContextPtr loc = this->getContextPtr();
    scaleRows.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<OtherType> rDiagonal( values, loc );
    WriteAccess<ValueType> wData( mData, loc );
    scaleRows[loc]( wData.get(), mNumRows, mNumColumns, rDiagonal.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::transposeImpl()
{
    SCAI_REGION( "Storage.Dense.transpose" )

    // Compute transpostion A^t of A via A^t = A^t * I, where * is implemented by LAPACK
    ContextPtr context = Context::getHostPtr();
    WriteAccess<ValueType> wData( mData, context );

    // transpose quadratic matrix
    // quadratic implementation is a special case of the rectangular one but this specific one might be faster
    if ( mNumColumns == mNumRows )
    {
        for ( IndexType i = 0; i < mNumColumns; ++i )
        {
            for ( IndexType j = i + 1; j < mNumColumns; ++j )
            {
                std::swap( wData[i + mNumColumns * j], wData[j + mNumColumns * i] );
            }
        }
    }
    else //tranpose rectangular matrix
    {
        for ( IndexType start = 0; start < mNumColumns * mNumRows; ++start )
        {
            IndexType next = start;
            IndexType i = 0;

            do
            {
                ++i;
                next = ( next % mNumRows ) * mNumColumns + next / mNumRows;
            }
            while ( next > start );

            if ( next >= start && i != 1 )
            {
                const ValueType tmp = wData[start];
                next = start;

                do
                {
                    i = ( next % mNumRows ) * mNumColumns + next / mNumRows;
                    wData[next] = ( i == start ) ? tmp :  wData[i];
                    next = i;
                }
                while ( next > start );
            }
        }
    }

    // swap dimensions of transposed matrix
    std::swap( mNumRows, mNumColumns );
};

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool DenseStorage<ValueType>::checkDiagonalProperty() const
{
    return mNumRows == mNumColumns;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t DenseStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( ValueType ) * mData.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType DenseStorage<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format DenseStorage<ValueType>::getFormat() const
{
    return Format::DENSE;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setIdentity( const IndexType size )
{
    if ( ( mNumRows == size ) && ( mNumColumns == size ) )
    {
        setZero();                            // no reallocation needed
    }
    else
    {
        allocate( size, size );               // initializes also with zero
    }

    setDiagonalImpl( ValueType( 1 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setZero()
{
    LAMAKernel<DenseKernelTrait::setValue<ValueType> > setValue;
    ContextPtr loc = this->getContextPtr();
    setValue.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<ValueType> data( mData, loc, mNumRows * mNumColumns );
    setValue[loc]( data.get(), mNumRows, mNumColumns, static_cast<ValueType>( 0 ), common::binary::COPY );
    SCAI_LOG_INFO( logger, *this << " has been set to zero" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorage<ValueType>::buildCSR(
    HArray<IndexType>& csrIA,
    HArray<IndexType>* csrJA,
    HArray<OtherValueType>* csrValues,
    const ContextPtr context ) const
{
    static LAMAKernel<DenseKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DenseKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;
    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;
    // check if context provides all implementations, otherwise go back to Host
    ContextPtr loc = context;
    getCSRValues.getSupportedContext( loc, getCSRSizes, sizes2offsets );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> denseValues( mData, loc );
    WriteOnlyAccess<IndexType> wIA( csrIA, loc, mNumRows + 1 );
    ValueType eps = this->mEpsilon;
    // Note: mDiagonalProperty == ( mNumRows == mNumColumns )
    getCSRSizes[loc]( wIA.get(), mDiagonalProperty, mNumRows, mNumColumns, denseValues.get(), eps );

    if ( csrJA == NULL || csrValues == NULL )
    {
        wIA.resize( mNumRows );
        return;
    }

    // build offset array, get number of non-zero values for size of ja, values
    IndexType numValues = sizes2offsets[loc]( wIA.get(), mNumRows );
    WriteOnlyAccess<IndexType> wJA( *csrJA, loc, numValues );
    WriteOnlyAccess<OtherValueType> wValues( *csrValues, loc, numValues );
    getCSRValues[loc]( wJA.get(), wValues.get(), wIA.get(), mDiagonalProperty, mNumRows, mNumColumns,
                       denseValues.get(), eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::getFirstColumnIndexes( hmemo::HArray<IndexType>& ) const
{
    COMMON_THROWEXCEPTION( "getFirstColumnIndexes not possible for DENSE format" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr context )
{
    static LAMAKernel<DenseKernelTrait::setCSRValues<ValueType, OtherValueType> > setCSRValues;
    static LAMAKernel<CSRKernelTrait::validOffsets> validOffsets;
    static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;

    // check if context provides all implementations, otherwise go back to Host

    ContextPtr loc = context;
    setCSRValues.getSupportedContext( loc, validOffsets, validIndexes );

    SCAI_LOG_INFO( logger,
                   "setCRSData for dense storage " << numRows << " x " << numColumns << ", nnz = " << numValues )
    mNumRows = numRows;
    mNumColumns = numColumns;
    std::unique_ptr<HArray<IndexType> > tmpOffsets;

    const HArray<IndexType>* offsets = &ia;

    if ( ia.size() == numRows )
    {
        tmpOffsets.reset( ia.copy() );
        IndexType total = HArrayUtils::scan1( *tmpOffsets, loc );
        SCAI_ASSERT_EQUAL( total, numValues, "sizes do not sum up correctly" )
        offsets = tmpOffsets.get();
    }
    else
    {
        SCAI_ASSERT_EQ_ERROR( ia.size(), numRows + 1, "size mismatch of csr IA array" )
    }

    {
        ReadAccess<IndexType> csrIA( *offsets, loc );
        ReadAccess<IndexType> csrJA( ja, loc );
        ReadAccess<OtherValueType> csrValues( values, loc );
        SCAI_CONTEXT_ACCESS( loc )

        if ( !validOffsets[loc]( csrIA.get(), numRows, numValues ) )
        {
            COMMON_THROWEXCEPTION( "invalid offset array" )
        }

        if ( !validIndexes[loc]( csrJA.get(), numValues, numColumns ) )
        {
            COMMON_THROWEXCEPTION( "invalid column indexes, #columns = " << numColumns )
        }

        WriteOnlyAccess<ValueType> data( mData, loc, mNumRows * mNumColumns );
        setCSRValues[loc]( data.get(), mNumRows, mNumColumns, csrIA.get(), csrJA.get(), csrValues.get() );
    }

    mDiagonalProperty = checkDiagonalProperty();
    // dense storage does not care about diagonal property, is always okay
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorage<ValueType>::setDIADataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const HArray<IndexType>& offsets,
    const HArray<OtherValueType>& values,
    const ContextPtr /* prefLoc */ )
{
    DIAStorage<OtherValueType> dia( numRows, numColumns, numDiagonals, offsets, values );
    assign( dia );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::invert( const MatrixStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "invert( " << other << ") to a dense storage" )

    if ( other.getFormat() == Format::DENSE )
    {
        const DenseStorage<ValueType>* otherDense = dynamic_cast<const DenseStorage<ValueType>*>( &other );
        SCAI_ASSERT_ERROR( otherDense, "Internal error: dynamic cast Dense" )
        invertDense( *otherDense );
    }
    else
    {
        SCAI_UNSUPPORTED( "invert (" << other << ") requires conversion to dense storage" )
        assign( other );
        invertDense( *this ); // is always done in place
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::invertDense( const DenseStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "invertDense: " << other )

    // invert is always done in place, so assign other to this

    if ( &other != this )
    {
        SCAI_LOG_INFO( logger, "invertDense: copy input matrix to this matrix" )
        assign( other );
    }

    int nRows = this->getNumRows();
    int nCols = this->getNumColumns();
    SCAI_ASSERT_EQUAL_ERROR( nRows, nCols )
    static LAMAKernel<blaskernel::BLASKernelTrait::getinv<ValueType> > getinv;
    ContextPtr loc = this->getContextPtr();
    getinv.getSupportedContext( loc );
    WriteAccess<ValueType> denseValues( this->getData(), loc );
    SCAI_CONTEXT_ACCESS( loc );
    getinv[loc]( nRows, denseValues.get(), nCols );
    SCAI_LOG_INFO( logger, "invertDense: this = " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y" << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    if ( alpha == common::constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::compute( result, beta, common::binary::MULT, y, this->getContextPtr() );
        return;
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )

    if ( beta != common::constants::ZERO )
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
    }

    if ( mNumRows == 0 )
    {
        result.clear();  // result will get size zero
        return;          // nothing to do
    }

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector, try on " << *mContext )

    // using BLAS2 interface requires result and y to be aliased

    if ( beta == common::constants::ZERO )
    {
        result.resize( mNumRows );
        utilskernel::HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, this->getContextPtr() );
    }
    else if ( &result != &y )
    {
        utilskernel::HArrayUtils::assign( result, y, this->getContextPtr() );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * A * x + beta * result

    if ( mNumColumns == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if ( beta == common::constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if ( beta == common::constants::ONE )
        {
            // no scaling required
        }
        else
        {
            HArrayUtils::compute( result, beta, common::binary::MULT, result, this->getContextPtr() );
        }
    }
    else
    {
        // mNumColums > 0, mnumRows > 0, so we avoid problems for gemv with m==0 or n==0
        static LAMAKernel<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
        ContextPtr loc = this->getContextPtr();
        gemv.getSupportedContext( loc );
        ReadAccess<ValueType> denseValues( mData, loc );
        ReadAccess<ValueType> rX( x, loc );
        int lda = mNumColumns; // stride for denseValues between rows
        WriteAccess<ValueType> wResult( result, loc );
        SCAI_CONTEXT_ACCESS( loc )
        // gemv:  result = alpha * this * x + beta * result
        gemv[loc]( CblasRowMajor, CblasNoTrans, mNumRows, mNumColumns, alpha, denseValues.get(), lda, rX.get(),
                   1, beta, wResult.get(), 1 );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger, "Computing z = " << alpha << " * x * A + " << beta << " * y"
                   << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )

    if ( beta != common::constants::ZERO )
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
    }

    if ( mNumColumns == 0 )
    {
        result.clear();  // result will get also zero size
        return;          // nothing more to do
    }

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector try on " << *mContext )

    // not used here: LAMA_INTERFACE_FN_T( normalGEVM, loc, DenseUtils, Mult, ValueType )

    // using BLAS2 interface requires result and y to be aliased

    if ( beta == common::constants::ZERO )
    {
        result.resize( mNumColumns );
        utilskernel::HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, this->getContextPtr() );
    }
    else if ( &result != &y )
    {
        SCAI_LOG_INFO( logger, "set result = y as y != result" )
        utilskernel::HArrayUtils::assign( result, y, this->getContextPtr() );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * x * A + beta * result

    if ( mNumRows == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if ( beta == common::constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if ( beta == common::constants::ONE )
        {
            // no scaling required
        }
        else
        {
            utilskernel::HArrayUtils::compute( result, result, common::binary::MULT, beta, this->getContextPtr() );
        }
    }
    else
    {
        // mNumColums > 0, mnumRows > 0, so we avoid problems for gevm with m==0 or n==0
        static LAMAKernel<blaskernel::BLASKernelTrait::gemv<ValueType> > gemv;
        ContextPtr loc = this->getContextPtr();
        gemv.getSupportedContext( loc );
        ReadAccess<ValueType> denseValues( mData, loc );
        ReadAccess<ValueType> rX( x, loc );
        int lda = mNumColumns; // stride for denseValues between rows
        WriteAccess<ValueType> wResult( result, loc );
        SCAI_CONTEXT_ACCESS( loc )
        // gemv:  result = alpha * this * x + beta * result
        gemv[loc]( CblasRowMajor, CblasTrans, mNumRows, mNumColumns, alpha, denseValues.get(), lda, rX.get(),
                   1, beta, wResult.get(), 1 );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    ReadAccess<ValueType> values( mData );
    stream << "DenseStorage " << mNumRows << " x " << mNumColumns << ", addr  = " << values.get() << endl;

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        stream << "Row " << i << " :";

        for ( IndexType j = 0; j < mNumColumns; j++ )
        {
            stream << " " << values[i * mNumColumns + j];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& c )
{
    SCAI_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B + " << beta << " * C" << ", with A = " << a << ", B = " << b << ", C = " << c )
    // a and b have to be Dense storages, otherwise create temporaries.
    const DenseStorage<ValueType>* denseA = NULL;
    const DenseStorage<ValueType>* denseB = NULL;
    const DenseStorage<ValueType>* denseC = NULL;
    // Define shared pointers in case we need temporaries
    shared_ptr<DenseStorage<ValueType> > tmpA;
    shared_ptr<DenseStorage<ValueType> > tmpB;
    shared_ptr<DenseStorage<ValueType> > tmpC;

    if ( a.getFormat() == Format::DENSE )
    {
        denseA = dynamic_cast<const DenseStorage<ValueType>*>( &a );
        SCAI_ASSERT_DEBUG( denseA, "could not cast to DenseStorage " << a )
    }
    else
    {
        SCAI_UNSUPPORTED( a << ": will be converted to Dense for matrix multiply" )
        tmpA = shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( a ) );
        denseA = tmpA.get();
    }

    if ( b.getFormat() == Format::DENSE )
    {
        denseB = dynamic_cast<const DenseStorage<ValueType>*>( &b );
        SCAI_ASSERT_DEBUG( denseB, "could not cast to DenseStorage " << b )
    }
    else
    {
        SCAI_UNSUPPORTED( b << ": will be converted to Dense for matrix multiply" )
        tmpB = shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( b ) );
        denseB = tmpB.get();
    }

    if ( c.getFormat() == Format::DENSE )
    {
        denseC = dynamic_cast<const DenseStorage<ValueType>*>( &c );
        SCAI_ASSERT_DEBUG( denseB, "could not cast to DenseStorage " << c )
    }
    else
    {
        SCAI_UNSUPPORTED( c << ": will ce converted to Dense for matrix multiply" )
        tmpC = shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( c ) );
        denseC = tmpC.get();
    }

    // now we have in any case all arguments as Dense Storage
    matrixTimesMatrixDense( alpha, *denseA, *denseB, beta, *denseC );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::matrixTimesMatrixDense(
    const ValueType alpha,
    const DenseStorage<ValueType>& a,
    const DenseStorage<ValueType>& b,
    const ValueType beta,
    const DenseStorage<ValueType>& c )
{
    // shape(a) = m x k,  shape(b) = k x n
    const DenseStorage<ValueType>* ptrA = &a;
    const DenseStorage<ValueType>* ptrB = &b;
    std::shared_ptr<DenseStorage<ValueType> > tmpA;
    std::shared_ptr<DenseStorage<ValueType> > tmpB;
    SCAI_LOG_INFO( logger,
                   "matrixTimesMatrixDense: " << alpha << " * a * b + " << beta << " * c, with a = " << a << ", b = " << b << ", c = " << c )

    if ( &a == this )
    {
        SCAI_LOG_INFO( logger, "temporary for A in A * B ( dense storages) needed" )
        tmpA = std::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( a ) );
        ptrA = tmpA.get();
    }

    if ( &b == this )
    {
        // avoid two temporaries
        if ( &a == &b )
        {
            SCAI_LOG_INFO( logger, "same temporary for A and B in A * B ( dense storages) needed" )
            ptrB = tmpA.get();
        }
        else
        {
            SCAI_LOG_INFO( logger, "temporary for B in A * B ( dense storages) needed" )
            tmpB = std::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( b ) );
            ptrB = tmpB.get();
        }
    }

    IndexType m = a.getNumRows();
    IndexType k = b.getNumRows();
    IndexType n = b.getNumColumns();

    SCAI_ASSERT_EQUAL_ERROR( k, a.getNumColumns() )
    mNumRows = m;
    mNumColumns = n;

    if ( beta == common::constants::ZERO )
    {
        // do not care at all about C as it might be any dummy, or aliased to result
        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );
        SCAI_LOG_INFO( logger, "init this result with 0, size = " << m * n )
        WriteOnlyAccess<ValueType> resAccess( getData(), loc, m * n );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( resAccess.get(), m * n, ValueType( 0 ), common::binary::COPY );
    }
    else if ( this != &c )
    {
        // force result = C
        SCAI_ASSERT_EQUAL_ERROR( m, c.getNumRows() )
        SCAI_ASSERT_EQUAL_ERROR( n, c.getNumColumns() )
        mNumRows = m;
        mNumColumns = n;
        static LAMAKernel<blaskernel::BLASKernelTrait::copy<ValueType> > copy;
        ContextPtr loc = this->getContextPtr();
        copy.getSupportedContext( loc );
        ReadAccess<ValueType> cAccess( c.getData(), loc );
        WriteOnlyAccess<ValueType> resAccess( getData(), loc, m * n );
        SCAI_LOG_TRACE( logger, "Copying: res = c " )
        SCAI_CONTEXT_ACCESS( loc )
        copy[loc]( n * m, cAccess.get(), 1, resAccess.get(), 1 );
    }
    else
    {
        SCAI_LOG_INFO( logger, "results is aliased with C as required for gemm, beta = " << beta )
    }

    // now C used for BLAS3 call of GEMM is this matrix
    int lda = a.getNumColumns();
    int ldb = b.getNumColumns();
    int ldc = mNumColumns;
    // Now  C = alpha * A * B + beta * C, can use GEMM of BLAS3
    SCAI_LOG_TRACE( logger,
                    "GEMM( CblasRowMajor, CblasNoTrans, CblasNoTrans, m:" << m << ", n:" << n << ", k:" << k << ", alpha:" << alpha << ", A , lda:" << lda << " ,B ,ldb:" << ldb << " ,beta:" << beta << " C ,ldc:" << ldc << " )" )

    if ( lda != 0 && n != 0 && m != 0 )
    {
        static LAMAKernel<blaskernel::BLASKernelTrait::gemm<ValueType> > gemm;
        ContextPtr loc = this->getContextPtr();
        gemm.getSupportedContext( loc );
        ReadAccess<ValueType> aAccess( ptrA->getData(), loc );
        ReadAccess<ValueType> bAccess( ptrB->getData(), loc );
        WriteAccess<ValueType> resAccess( getData(), loc );
        SCAI_CONTEXT_ACCESS( loc )
        gemm[loc]( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, aAccess.get(), lda, bAccess.get(), ldb, beta,
                   resAccess.get(), ldc );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::~DenseStorage()
{
    SCAI_LOG_DEBUG( logger, "~DenseStorage for matrix " << mNumRows << " x " << mNumColumns )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseStorage<ValueType>::l1Norm() const
{
    IndexType n = mNumRows * mNumColumns;

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    ContextPtr loc = this->getContextPtr();
    asum.getSupportedContext( loc );
    ReadAccess<ValueType> data( mData, loc );
    SCAI_CONTEXT_ACCESS( loc );
    IndexType inc = 1;
    return asum[loc]( n, data.get(), inc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseStorage<ValueType>::l2Norm() const
{
    IndexType n = mNumRows * mNumColumns;

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = this->getContextPtr();
    dot.getSupportedContext( loc );
    ReadAccess<ValueType> data( mData, loc );
    SCAI_CONTEXT_ACCESS( loc );
    return common::Math::sqrt( dot[loc]( n, data.get(), 1, data.get(), 1 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseStorage<ValueType>::maxNorm() const
{
    IndexType n = mNumRows * mNumColumns;

    if ( n == 0 )
    {
        return ValueType( 0 );
    }

    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;
    ContextPtr loc = this->getContextPtr();
    reduce.getSupportedContext( loc );
    ReadAccess<ValueType> read1( mData, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ValueType zero   = 0;
    NormType<ValueType> maxval = reduce[loc]( read1.get(), n, zero, common::binary::ABS_MAX );
    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
NormType<ValueType> DenseStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows, other.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, other.getNumColumns() )
    std::shared_ptr<DenseStorage<ValueType> > tmpOtherDense;
    const DenseStorage<ValueType>* otherDense;

    if ( other.getValueType() == getValueType() && ( other.getFormat() == Format::DENSE ) )
    {
        otherDense = dynamic_cast<const DenseStorage<ValueType>*>( &other );
        SCAI_ASSERT_ERROR( otherDense, other << ": could not cast to " << typeName() )
    }
    else
    {
        SCAI_UNSUPPORTED( other << ": converted to " << typeName() << " for maxDiffNorm" )
        tmpOtherDense.reset( new DenseStorage<ValueType>( other ) );
        otherDense = tmpOtherDense.get();
    }

    return maxDiffNormImpl( *otherDense );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorage<ValueType>::maxDiffNormImpl( const DenseStorage<ValueType>& other ) const
{
    // no more checks needed here

    IndexType n = mNumRows * mNumColumns;

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
    }

    return HArrayUtils::absMaxDiffVal( mData, other.mData, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorage<ValueType>::assignDenseStorageImpl( const DenseStorage<OtherValueType>& other )
{
    // actualize member variables of base class
    _MatrixStorage::_assign( other ); // copy sizes, flags
    LAMAKernel<DenseKernelTrait::set<ValueType, OtherValueType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    {
        SCAI_CONTEXT_ACCESS( loc )
        WriteOnlyAccess<ValueType> data( mData, loc, mNumRows * mNumColumns );
        ReadAccess<OtherValueType> otherData( other.getData(), loc );
        set[loc]( data.get(), mNumRows, mNumColumns, otherData.get(), common::binary::COPY );
    }
    SCAI_LOG_INFO( logger, *this << ": assigned dense storage " << other )
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // Here more efficient solutions can be used instead of build/set CSR data
    if ( &other == this )
    {
        SCAI_LOG_INFO( logger, "self assign, is skipped" )
        return;
    }

    if ( other.getFormat() == Format::DENSE )
    {
        // more efficient solution for assigment of dense storage
        if ( mepr::DenseStorageWrapper<ValueType, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( *this, other ) )
        {
            return;
        }
    }

    SCAI_LOG_INFO( logger, *this << ": (dense) assign " << other )

    // In all other cases we use the fallback routine

    MatrixStorage<ValueType>::assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::allocate( IndexType numRows, IndexType numCols )
{
    SCAI_LOG_INFO( logger, "allocate Dense sparse matrix of size " << numRows << " x " << numCols )
    mNumRows = numRows;
    mNumColumns = numCols;
    setZero();
    SCAI_LOG_DEBUG( logger, *this << " allocated, #values = " << mNumRows * mNumColumns << ", initialized with 0" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << getTypeName() << "( rows = " << mNumRows << ", cols = " << mNumColumns << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    const IndexType pos = i * mNumColumns + j;

    return utilskernel::HArrayUtils::getVal<ValueType>( mData, pos );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::setValue( const IndexType i,
                                        const IndexType j,
                                        const ValueType val,
                                        const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    const IndexType pos = i * mNumColumns + j;

    utilskernel::HArrayUtils::setValImpl( mData, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mData.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::wait() const
{
    mData.wait();
}

/* --------------------------------------------------------------------------- */
/*  Constructors for DenseStorage                                              */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage() :

    CRTPMatrixStorage<DenseStorage<ValueType>, ValueType>( 0, 0 )
{
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const IndexType numRows, const IndexType numColumns ) :

    CRTPMatrixStorage<DenseStorage<ValueType>, ValueType>( numRows, numColumns )

{
    mDiagonalProperty = checkDiagonalProperty();
    this->setZero();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage(
    const HArray<ValueType>& data,
    const IndexType numRows,
    const IndexType numColumns ) :

    CRTPMatrixStorage<DenseStorage<ValueType>, ValueType>( numRows, numColumns )

{
    SCAI_ASSERT_EQUAL_ERROR( data.size(), numRows * numColumns )
    mData = data; // is a flat copy of the array
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const DenseStorage<ValueType>& other ) :

    CRTPMatrixStorage<DenseStorage<ValueType>, ValueType>()

{
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const _MatrixStorage& other ) :

    CRTPMatrixStorage<DenseStorage<ValueType>, ValueType>()

{
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const _MatrixStorage& other, const hmemo::ContextPtr loc ) :

    CRTPMatrixStorage<DenseStorage<ValueType>, ValueType>()

{
    _MatrixStorage::setContextPtr( loc );
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>& DenseStorage<ValueType>::operator=( const _MatrixStorage& other )
{
    assign( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>& DenseStorage<ValueType>::operator=( const DenseStorage<ValueType>& other )
{
    assign( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::swap( _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( getFormat(), other.getFormat(), "swap only for same storage format" )
    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(), "swap only for same value type" )

    // only in debug mode use the more expensive dynamic cast for verification

    SCAI_ASSERT_DEBUG( dynamic_cast<DenseStorage<ValueType>* >( &other ), "illegal storage to swap" )

    swapImpl( reinterpret_cast<DenseStorage<ValueType>& >( other ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::swapImpl( DenseStorage<ValueType>& other )
{
    MatrixStorage<ValueType>::swapMS( other ); // swap member variable of base class

    mData.swap( other.mData ); // swap data

    // other properties remain unchanged, also the context
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorage<ValueType>::swap( HArray<ValueType>& other, const IndexType numRows, const IndexType numColumns )
{
    SCAI_ASSERT_EQ_ERROR( other.size(), numRows * numColumns,
                          "array of dense storage has not size " << numRows << " x " << numColumns )

    mNumRows    = numRows;
    mNumColumns = numColumns;

    mData.swap( other );

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const char* DenseStorage<ValueType>::getTypeName() const
{
    return typeName();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>*
DenseStorage<ValueType>::copy() const
{
    return new DenseStorage<ValueType>( *this );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType DenseStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DENSE, common::getScalarType<ValueType>() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* DenseStorage<ValueType>::create()
{
    return new DenseStorage<ValueType>();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
DenseStorage<ValueType>* DenseStorage<ValueType>::newMatrixStorage() const
{
    std::unique_ptr<DenseStorage<ValueType> > storage( new DenseStorage<ValueType>() );
    storage->setContextPtr( this->getContextPtr() );
    return storage.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
std::string DenseStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "DenseStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const char* DenseStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseStorage, SCAI_NUMERIC_TYPES_HOST )

#define DENSE_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                             \
    template void DenseStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,        \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                            \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                             \
    template void DenseStorage<ValueType>::setDIADataImpl( const IndexType, const IndexType, const IndexType,        \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );            \
    template void DenseStorage<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;      \
    template void DenseStorage<ValueType>::setRowImpl( const hmemo::HArray<OtherValueType>&, const IndexType,        \
            const common::binary::BinaryOp );                  \
    template void DenseStorage<ValueType>::setColumnImpl( const hmemo::HArray<OtherValueType>&, const IndexType,     \
            const common::binary::BinaryOp );               \
    template void DenseStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                  \
    template void DenseStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                  \
    template void DenseStorage<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                        \
    template void DenseStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,           \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;  \
     
#define DENSE_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, DENSE_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( DENSE_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef DENSE_STORAGE_INST_LVL2
#undef DENSE_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
