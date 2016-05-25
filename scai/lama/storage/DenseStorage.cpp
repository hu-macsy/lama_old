/**
 * @file DenseStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Instantiation for template class DenseStorage.
 * @author Thomas Brandes, Michael Drost
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/DenseStorage.hpp>

// local libraries
#include <scai/lama/mepr/DenseStorageViewWrapper.hpp>

// internal scai libraries
#include <scai/sparsekernel/DenseKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/hmemo/ContextAccess.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::hmemo;

using std::abs;
// so we can use abs for float and double and own abs for Complex

namespace scai
{

using common::shared_ptr;
using common::TypeTraits;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using utilskernel::UtilKernelTrait;

using sparsekernel::DenseKernelTrait;
using sparsekernel::CSRKernelTrait;


namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseStorageView<ValueType>::logger,
                              "MatrixStorage.DenseStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorageView<ValueType>::DenseStorageView(
    HArray<ValueType>& data,
    const IndexType numRows,
    const IndexType numColumns,
    bool initializedData )

    : CRTPMatrixStorage<DenseStorageView<ValueType>, ValueType>( numRows, numColumns ), mData( data )
{
    if ( initializedData )
    {
        SCAI_ASSERT_EQUAL_ERROR( data.size(), numRows * numColumns )
    }

    SCAI_LOG_DEBUG( logger, "created dense storage with referenced data" )

    // we can always access diagonals directly

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
HArray<ValueType>& DenseStorageView<ValueType>::getData()
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& DenseStorageView<ValueType>::getData() const
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DenseStorageView<ValueType>::getNumValues() const
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
void DenseStorageView<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void DenseStorageView<ValueType>::check( const char* /* msg */ ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows * mNumColumns, mData.size() )
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::setDiagonalImpl( const ValueType value )
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
template<typename OtherType>
void DenseStorageView<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    static LAMAKernel<DenseKernelTrait::getRow<OtherType, ValueType> > getRow;

    ContextPtr loc = this->getContextPtr();
    getRow.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<OtherType> wRow( row, loc, mNumColumns );
    ReadAccess<ValueType> rData( mData, loc );

    getRow[loc]( wRow.get(), rData.get(), i, mNumRows, mNumColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::getDiagonalImpl( HArray<OtherType>& diagonal ) const
{
    IndexType numDiagonalValues = std::min( mNumColumns, mNumRows );

    static LAMAKernel<DenseKernelTrait::getDiagonal<OtherType, ValueType> > getDiagonal;

    ContextPtr loc = this->getContextPtr();
    getDiagonal.getSupportedContext( loc );

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalValues );
    ReadAccess<ValueType> rData( mData, loc );

    SCAI_CONTEXT_ACCESS( loc )

    getDiagonal[loc]( wDiagonal.get(), numDiagonalValues, rData.get(), mNumRows, mNumColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::setDiagonalImpl( const HArray<OtherType>& diagonal )
{
    IndexType numDiagonalValues = std::min( mNumColumns, mNumRows );

    static LAMAKernel<DenseKernelTrait::setDiagonal<ValueType, OtherType> > setDiagonal;

    ContextPtr loc = this->getContextPtr();
    setDiagonal.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wData( mData, loc );

    setDiagonal[loc]( wData.get(), mNumRows, mNumColumns, rDiagonal.get(), numDiagonalValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::scaleImpl( const ValueType value )
{
   
    // not used here: HArrayUtils::scale( mData, value, this->getContextPtr() )
    // reasoning:     workload distribution would not fit to distribution of rows

    static LAMAKernel<DenseKernelTrait::scaleValue<ValueType> > scaleValue;

    ContextPtr loc = this->getContextPtr();

    scaleValue.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wData( mData, loc );

    scaleValue[loc]( wData.get(), mNumRows, mNumColumns, value );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::conj()
{
    HArrayUtils::conj( mData, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::scaleImpl( const HArray<OtherType>& values )
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
void DenseStorageView<ValueType>::transposeImpl(){
    // Compute transpostion A^t of A via A^t = A^t * I, where * is implemented by LAPACK
    ContextPtr context = Context::getHostPtr();
    WriteAccess<ValueType> wData( mData, context );

    // transpose quadratic matrix
    // quadratic implementation is a special case of the rectangular one but this specific one might be faster
    if(mNumColumns == mNumRows)
    {
        for(IndexType i=0;i<mNumColumns;++i)
        {
            for(IndexType j=i+1; j<mNumColumns;++j)
                std::swap(wData[i+mNumColumns*j], wData[j+mNumColumns*i]);
        }  
    }
    else //tranpose rectangular matrix
    {
        for(IndexType start=0; start< mNumColumns*mNumRows;++start)
        {
            IndexType next = start;
            IndexType i = 0;
            do{
                ++i;
                next= (next % mNumRows) * mNumColumns + next / mNumRows;
            } while(next > start);

            if(next >= start && i != 1)
            {
                const ValueType tmp= wData[start];
                next = start;
                do{
                    i = (next % mNumRows) * mNumColumns + next / mNumRows;
                    wData[next] = (i == start ) ? tmp :  wData[i];
                    next = i;
                } while(next > start);
            }
        }
    }

    // swap dimensions of transposed matrix

    std::swap( mNumRows, mNumColumns );
};

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool DenseStorageView<ValueType>::checkDiagonalProperty() const
{
    return mNumRows == mNumColumns;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t DenseStorageView<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( ValueType ) * mData.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
common::scalar::ScalarType DenseStorageView<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format::MatrixStorageFormat DenseStorageView<ValueType>::getFormat() const
{
    return Format::DENSE;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::setIdentity( const IndexType size )
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
void DenseStorageView<ValueType>::setZero()
{
    LAMAKernel<DenseKernelTrait::setValue<ValueType> > setValue;

    ContextPtr loc = this->getContextPtr();
    setValue.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<ValueType> data( mData, loc, mNumRows * mNumColumns );

    setValue[loc]( data.get(), mNumRows, mNumColumns, static_cast<ValueType>( 0 ) );

    SCAI_LOG_INFO( logger, *this << " has been set to zero" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorageView<ValueType>::buildCSR(
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
template<typename OtherValueType>
void DenseStorageView<ValueType>::setCSRDataImpl(
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
    setCSRValues.getSupportedContext( loc );

    SCAI_LOG_INFO( logger,
                   "setCRSData for dense storage " << numRows << " x " << numColumns << ", nnz = " << numValues )

    mNumRows = numRows;
    mNumColumns = numColumns;

    common::unique_ptr<HArray<IndexType> > tmpOffsets;

    const HArray<IndexType>* offsets = &ia;

    if ( ia.size() == numRows )
    {
        tmpOffsets.reset( ia.copy() );
        IndexType total = HArrayUtils::scan( *tmpOffsets, loc );
        SCAI_ASSERT_EQUAL( total, numValues, "sizes do not sum up correctly" )
        offsets = tmpOffsets.get();
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
void DenseStorageView<ValueType>::invert( const MatrixStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "invert( " << other << ") to a dense storage" )

    if ( other.getFormat() == Format::DENSE )
    {
        const DenseStorageView<ValueType>* otherDense = dynamic_cast<const DenseStorageView<ValueType>*>( &other );

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
void DenseStorageView<ValueType>::invertDense( const DenseStorageView<ValueType>& other )
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
void DenseStorageView<ValueType>::matrixTimesVector(
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

        HArrayUtils::assignScaled( result, beta, y, this->getContextPtr() );
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

    if ( beta == scai::common::constants::ZERO )
    {
        SCAI_LOG_INFO( logger, "set result = 0 as y != result and beta = 0" )

        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), mNumRows, ValueType( 0 ), utilskernel::reduction::COPY );
    }
    else if ( &result != &y )
    {
        SCAI_LOG_INFO( logger, "set result = y as y != result" )

        static LAMAKernel<blaskernel::BLASKernelTrait::copy<ValueType> > copy;

        ContextPtr loc = this->getContextPtr();
        copy.getSupportedContext( loc );

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

        // by setting result = y we get the correct results

        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        copy[loc]( mNumRows, rY.get(), 1, wResult.get(), 1 );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * A * x + beta * result

    if ( mNumColumns == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if ( beta == scai::common::constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if ( beta == scai::common::constants::ONE )
        {
            // no scaling required
        }
        else
        {
            static LAMAKernel<blaskernel::BLASKernelTrait::scal<ValueType> > scal;

            ContextPtr loc = this->getContextPtr();
            scal.getSupportedContext( loc );

            WriteAccess<ValueType> wResult( result, loc );

            SCAI_CONTEXT_ACCESS( loc )

            scal[loc]( mNumRows, beta, wResult.get(), 1 );
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
void DenseStorageView<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y" << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )

    if ( beta != common::constants::ZERO )
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
    }

    if ( mNumRows == 0 )
    {
        return; // nothing to do
    }

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector try on " << *mContext )

    // not used here: LAMA_INTERFACE_FN_T( normalGEVM, loc, DenseUtils, Mult, ValueType )

    // using BLAS2 interface requires result and y to be aliased

    if ( beta == scai::common::constants::ZERO )
    {
        SCAI_LOG_INFO( logger, "set result = 0 as y != result and beta = 0" )

        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), mNumColumns, ValueType( 0 ), utilskernel::reduction::COPY );
    }
    else if ( &result != &y )
    {
        SCAI_LOG_INFO( logger, "set result = y as y != result" )

        static LAMAKernel<blaskernel::BLASKernelTrait::copy<ValueType> > copy;

        ContextPtr loc = this->getContextPtr();
        copy.getSupportedContext( loc );

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

        // by setting result = y we get the correct results

        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        copy[loc]( mNumColumns, rY.get(), 1, wResult.get(), 1 );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * A * x + beta * result

    if ( mNumColumns == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if ( beta == scai::common::constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if ( beta == scai::common::constants::ONE )
        {
            // no scaling required
        }
        else
        {
            static LAMAKernel<blaskernel::BLASKernelTrait::scal<ValueType> > scal;

            ContextPtr loc = this->getContextPtr();
            scal.getSupportedContext( loc );

            WriteAccess<ValueType> wResult( result, loc );

            SCAI_CONTEXT_ACCESS( loc )

            scal[loc]( mNumColumns, beta, wResult.get(), 1 );
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
void DenseStorageView<ValueType>::print( std::ostream& stream ) const
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
void DenseStorageView<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& c )
{
    SCAI_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B + " << beta << " * C" << ", with A = " << a << ", B = " << b << ", C = " << c )

    // a and b have to be Dense storages, otherwise create temporaries.

    const DenseStorageView<ValueType>* denseA = NULL;
    const DenseStorageView<ValueType>* denseB = NULL;
    const DenseStorageView<ValueType>* denseC = NULL;

    // Define shared pointers in case we need temporaries

    shared_ptr<DenseStorage<ValueType> > tmpA;
    shared_ptr<DenseStorage<ValueType> > tmpB;
    shared_ptr<DenseStorage<ValueType> > tmpC;

    if ( a.getFormat() == Format::DENSE )
    {
        denseA = dynamic_cast<const DenseStorageView<ValueType>*>( &a );
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
        denseB = dynamic_cast<const DenseStorageView<ValueType>*>( &b );
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
        denseC = dynamic_cast<const DenseStorageView<ValueType>*>( &c );
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
void DenseStorageView<ValueType>::matrixTimesMatrixDense(
    const ValueType alpha,
    const DenseStorageView<ValueType>& a,
    const DenseStorageView<ValueType>& b,
    const ValueType beta,
    const DenseStorageView<ValueType>& c )
{
    // shape(a) = m x k,  shape(b) = k x n

    const DenseStorageView<ValueType>* ptrA = &a;
    const DenseStorageView<ValueType>* ptrB = &b;

    common::shared_ptr<DenseStorage<ValueType> > tmpA;
    common::shared_ptr<DenseStorage<ValueType> > tmpB;

    SCAI_LOG_INFO( logger,
                   "matrixTimesMatrixDense: " << alpha << " * a * b + " << beta << " * c, with a = " << a << ", b = " << b << ", c = " << c )

    if ( &a == this )
    {
        SCAI_LOG_INFO( logger, "temporary for A in A * B ( dense storages) needed" )
        tmpA = common::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( a ) );
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
            tmpB = common::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( b ) );
            ptrB = tmpB.get();
        }
    }

    int m = a.getNumRows();
    int k = b.getNumRows();
    int n = b.getNumColumns();

    SCAI_ASSERT_EQUAL_ERROR( k, a.getNumColumns() )

    mNumRows = m;
    mNumColumns = n;

    if ( beta == scai::common::constants::ZERO )
    {
        // do not care at all about C as it might be any dummy, or aliased to result

        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );

        SCAI_LOG_INFO( logger, "init this result with 0, size = " << m * n )
        WriteOnlyAccess<ValueType> resAccess( getData(), loc, m * n );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( resAccess.get(), m * n, ValueType( 0 ), utilskernel::reduction::COPY );
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
DenseStorageView<ValueType>::~DenseStorageView()
{
    SCAI_LOG_DEBUG( logger, "~DenseStorage for matrix " << mNumRows << " x " << mNumColumns )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorageView<ValueType>::l1Norm() const
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
ValueType DenseStorageView<ValueType>::l2Norm() const
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
ValueType DenseStorageView<ValueType>::maxNorm() const
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

    ValueType maxval = reduce[loc]( read1.get(), n, utilskernel::reduction::ABS_MAX );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorageView<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows, other.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, other.getNumColumns() )

    common::shared_ptr<DenseStorage<ValueType> > tmpOtherDense;

    const DenseStorageView<ValueType>* otherDense;

    if ( other.getValueType() == getValueType() && ( other.getFormat() == Format::DENSE ) )
    {
        otherDense = dynamic_cast<const DenseStorageView<ValueType>*>( &other );
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
ValueType DenseStorageView<ValueType>::maxDiffNormImpl( const DenseStorageView<ValueType>& other ) const
{
    // no more checks needed here

    IndexType n = mNumRows * mNumColumns;

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
    }

    static LAMAKernel<UtilKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;

    ContextPtr loc = this->getContextPtr();
    absMaxDiffVal.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> read1( mData, loc );
    ReadAccess<ValueType> read2( other.mData, loc );

    ValueType maxval = absMaxDiffVal[loc]( read1.get(), read2.get(), n );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorageView<ValueType>::assignDenseStorageImpl( const DenseStorageView<OtherValueType>& other )
{
    // actualize member variables of base class

    _MatrixStorage::_assign( other ); // copy sizes, flags

    LAMAKernel<DenseKernelTrait::copyDenseValues<ValueType, OtherValueType> > copyDenseValues;

    ContextPtr loc = this->getContextPtr();
    copyDenseValues.getSupportedContext( loc );

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<ValueType> data( mData, loc, mNumRows * mNumColumns );
        ReadAccess<OtherValueType> otherData( other.getData(), loc );

        copyDenseValues[loc]( data.get(), mNumRows, mNumColumns, otherData.get() );
    }

    SCAI_LOG_INFO( logger, *this << ": assigned dense storage " << other )

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::assign( const _MatrixStorage& other )
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

        if( mepr::DenseStorageViewWrapper<ValueType, SCAI_ARITHMETIC_HOST_LIST>::assignImpl( *this, other ) )
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
void DenseStorageView<ValueType>::allocate( IndexType numRows, IndexType numCols )
{
    SCAI_LOG_INFO( logger, "allocate Dense sparse matrix of size " << numRows << " x " << numCols )

    mNumRows = numRows;
    mNumColumns = numCols;

    setZero();

    SCAI_LOG_DEBUG( logger, *this << " allocated, #values = " << mNumRows * mNumColumns << ", initialized with 0" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << getTypeName() << "( rows = " << mNumRows << ", cols = " << mNumColumns << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorageView<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    const ReadAccess<ValueType> data( mData );
    SCAI_LOG_TRACE( logger,
                    "get value (" << i << ", " << j << ")--> " << i* mNumColumns + j << "= " << data[i * mNumColumns + j] )
    return data[i * mNumColumns + j];
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::prefetch( const ContextPtr location ) const
{
    mData.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::wait() const
{
    mData.wait();
}

/* --------------------------------------------------------------------------- */
/*  Constructors for DenseStorage                                              */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage()

    : DenseStorageView<ValueType>( mDataArray, 0, 0, false )
{
// mDataArray is now correctly initialized
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const IndexType numRows, const IndexType numColumns )

    : DenseStorageView<ValueType>( mDataArray, 0, 0, false )
{
    // now resize the array and fill it with zero

    mNumRows = numRows;
    mNumColumns = numColumns;

    this->setZero();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage(
    const HArray<ValueType>& data,
    const IndexType numRows,
    const IndexType numColumns )

    : DenseStorageView<ValueType>( mDataArray, numRows, numColumns, false )

{
// here we can check for correct size

    SCAI_ASSERT_EQUAL_ERROR( data.size(), numRows * numColumns )

    mDataArray = data; // is a flat copy of the array
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const DenseStorage<ValueType>& other )

    : DenseStorageView<ValueType>( mDataArray, 0, 0, false )

{
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const _MatrixStorage& other )

    : DenseStorageView<ValueType>( mDataArray, 0, 0, false )

{
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage( const _MatrixStorage& other, const hmemo::ContextPtr loc )

    : DenseStorageView<ValueType>( mDataArray, 0, 0, false )

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
void DenseStorageView<ValueType>::swap( DenseStorageView<ValueType>& other )
{
    MatrixStorage<ValueType>::swap( other ); // swap dimensions
    mData.swap( other.mData ); // swap data
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::~DenseStorage()
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const char* DenseStorage<ValueType>::getTypeName() const
{
    return typeName();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorageView<ValueType>*
DenseStorageView<ValueType>::copy() const
{
    return new DenseStorage<ValueType>( *this );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* DenseStorageView<ValueType>::create()
{
    COMMON_THROWEXCEPTION( "creation currently not possible")
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType DenseStorageView<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DENSE, common::getScalarType<ValueType>() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
DenseStorageView<ValueType>* DenseStorageView<ValueType>::newMatrixStorage() const
{
    COMMON_THROWEXCEPTION( "DenseStorageView<ValueType>::newMatrixStorage() not implemented yet" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* DenseStorage<ValueType>::create()
{
    return new DenseStorage<ValueType>();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType DenseStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DENSE, common::getScalarType<ValueType>() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
DenseStorage<ValueType>* DenseStorage<ValueType>::newMatrixStorage() const
{
   common::unique_ptr<DenseStorage<ValueType> > storage( new DenseStorage<ValueType>() ); 
   storage->setContextPtr( this->getContextPtr() );
   return storage.release();
}

template<typename ValueType>
std::string DenseStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string("DenseStorage<") << common::getScalarType<ValueType>() << std::string(">");
    return s.str();
}

template<typename ValueType>
const char* DenseStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
std::string DenseStorageView<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string("DenseStorageView<") << common::getScalarType<ValueType>() << std::string(">");
    return s.str();
}

template<typename ValueType>
const char* DenseStorageView<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DenseStorage, SCAI_ARITHMETIC_HOST )
SCAI_COMMON_INST_CLASS( DenseStorageView, SCAI_ARITHMETIC_HOST )

#define DENSE_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                  \
     template void DenseStorageView<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,                \
                                                          const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&, \
                                                          const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );  \
     template void DenseStorageView<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;              \
     template void DenseStorageView<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                          \
     template void DenseStorageView<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                          \
     template void DenseStorageView<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                                \
     template void DenseStorageView<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
                                                    hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;  \

#define DENSE_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, DENSE_STORAGE_INST_LVL2, SCAI_ARITHMETIC_HOST )

SCAI_COMMON_LOOP( DENSE_STORAGE_INST_LVL1, SCAI_ARITHMETIC_HOST )

#undef DENSE_STORAGE_INST_LVL2
#undef DENSE_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
