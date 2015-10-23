/**
 * @file DenseStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @brief Instantiation for template class DenseStorage.
 * @author Thomas Brandes, Michael Drost
 * @date 04.06.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/storage/DenseStorage.hpp>

// local library
#include <scai/lama/BLASInterface.hpp>
#include <scai/lama/UtilsInterface.hpp>
#include <scai/lama/LAMAKernel.hpp>

#include <scai/lama/openmp/OpenMPDenseUtils.hpp>
#include <scai/lama/openmp/OpenMPCSRUtils.hpp>
#include <scai/lama/openmp/OpenMPUtils.hpp>

// internal scai libraries
#include <scai/hmemo/ContextAccess.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/exception/UnsupportedException.hpp>

using namespace scai::hmemo;

using std::abs;
// so we can use abs for float and double and own abs for Complex

namespace scai
{

using common::shared_ptr;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DenseStorageView<ValueType>::logger,
                              "MatrixStorage.DenseStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorageView<ValueType>::DenseStorageView(
    LAMAArray<ValueType>& data,
    const IndexType numRows,
    const IndexType numColumns,
    bool initializedData )

    : CRTPMatrixStorage<DenseStorageView<ValueType>,ValueType>( numRows, numColumns ), mData( data )
{
    if( initializedData )
    {
        SCAI_ASSERT_EQUAL_ERROR( data.size(), numRows * numColumns )
    }

    SCAI_LOG_DEBUG( logger, "created dense storage with referenced data" )

    // we can always access diagonals directly

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
LAMAArray<ValueType>& DenseStorageView<ValueType>::getData()
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<ValueType>& DenseStorageView<ValueType>::getData() const
{
    return mData;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DenseStorageView<ValueType>::getNumValues() const
{
    ReadAccess<ValueType> values( mData );

    IndexType count = 0;

    for( IndexType i = 0; i < mNumRows; ++i )
    {
        for( IndexType j = 0; j < mNumColumns; ++j )
        {
            if( abs( values[i * mNumColumns + j] ) > MatrixStorage<ValueType>::mEpsilon )
            {
                count++;
            }
        }
    }

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
void DenseStorageView<ValueType>::check( const char* /* msg */) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows * mNumColumns, mData.size() )
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::setDiagonalImpl( const Scalar value )
{
    WriteAccess<ValueType> wData( mData ); // use existing data

    OpenMPDenseUtils::setDiagonalValue( wData.get(), mNumRows, mNumColumns, value.getValue<ValueType>() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    WriteOnlyAccess<OtherType> wRow( row, mNumColumns );
    ReadAccess<ValueType> rData( mData );

    #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

    for( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = static_cast<OtherType>( rData[i * mNumColumns + j] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const
{
    IndexType numDiagonalValues = std::min( mNumColumns, mNumRows );

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, numDiagonalValues );

    ReadAccess<ValueType> rData( mData );

    OpenMPDenseUtils::getDiagonal( wDiagonal.get(), numDiagonalValues, rData.get(), mNumRows, mNumColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::setDiagonalImpl( const LAMAArray<OtherType>& diagonal )
{
    IndexType numDiagonalValues = std::min( mNumColumns, mNumRows );
    ReadAccess<OtherType> rDiagonal( diagonal );
    WriteAccess<ValueType> wData( mData );

    OpenMPDenseUtils::setDiagonal( wData.get(), mNumRows, mNumColumns, rDiagonal.get(), numDiagonalValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::scaleImpl( const Scalar values )
{
    WriteAccess<ValueType> wData( mData );

    const ValueType val = values.getValue<ValueType>();

    OpenMPDenseUtils::scaleValue( wData.get(), mNumRows, mNumColumns, val );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DenseStorageView<ValueType>::scaleImpl( const LAMAArray<OtherType>& values )
{
    ReadAccess<OtherType> rDiagonal( values );
    WriteAccess<ValueType> wData( mData );

    const IndexType columns = mNumColumns;

    #pragma omp parallel for schedule (SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < mNumRows; ++i )
    {
        const ValueType tmp = static_cast<ValueType>( rDiagonal[i] );
        const IndexType offset = i * columns;

        for( IndexType j = 0; j < columns; ++j )
        {
            wData[offset + j] *= tmp;
        }
    }
}

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
    memoryUsage += sizeof(ValueType) * mData.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType DenseStorageView<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat DenseStorageView<ValueType>::getFormat() const
{
    return Format::DENSE;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::setIdentity( const IndexType size )
{
    allocate( size, size ); // does not initialize
    setIdentity();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::setIdentity()
{
    WriteOnlyAccess<ValueType> data( mData, mNumRows * mNumColumns );

    //  do not use something like scale, data can be completely undefined

    OpenMPDenseUtils::scaleValue( data.get(), mNumRows, mNumColumns, static_cast<ValueType>(0.0) );
    OpenMPDenseUtils::setDiagonalValue( data.get(), mNumRows, mNumColumns, static_cast<ValueType>(1.0) );

    SCAI_LOG_INFO( logger, *this << " has been set to identity" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::setZero()
{
    WriteOnlyAccess<ValueType> data( mData, mNumRows * mNumColumns );

    const ValueType zero = static_cast<ValueType>(0.0);

    //  do not use something like scale, data can be completely undefined

    OpenMPDenseUtils::scaleValue( data.get(), mNumRows, mNumColumns, zero );

    SCAI_LOG_INFO( logger, *this << " has been set to zero" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorageView<ValueType>::buildCSR(
    LAMAArray<IndexType>& csrIA,
    LAMAArray<IndexType>* csrJA,
    LAMAArray<OtherValueType>* csrValues,
    const ContextPtr /* loc */) const
{
    // TODO all done on host, so loc is unused

    ReadAccess<ValueType> denseValues( mData );

    WriteOnlyAccess<IndexType> wIA( csrIA, mNumRows + 1 );

    ValueType eps = this->mEpsilon;

    // Note: mDiagonalProperty == ( mNumRows == mNumColumns )

    OpenMPDenseUtils::getCSRSizes( wIA.get(), mDiagonalProperty, mNumRows, mNumColumns, denseValues.get(), eps );

    if( csrJA == NULL || csrValues == NULL )
    {
        wIA.resize( mNumRows );
        return;
    }

    // build offset array, get number of non-zero values for size of ja, values

    IndexType numValues = OpenMPCSRUtils::sizes2offsets( wIA.get(), mNumRows );

    WriteOnlyAccess<IndexType> wJA( *csrJA, numValues );
    WriteOnlyAccess<OtherValueType> wValues( *csrValues, numValues );

    OpenMPDenseUtils::getCSRValues( wJA.get(), wValues.get(), wIA.get(), mDiagonalProperty, mNumRows, mNumColumns,
                                    denseValues.get(), eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DenseStorageView<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    const ContextPtr /* loc */)
{
    // not yet suppored on other devices

    ContextPtr loc = Context::getContextPtr( context::Host );

    SCAI_LOG_INFO( logger,
                   "setCRSData for dense storage " << numRows << " x " << numColumns << ", nnz = " << numValues )

    mNumRows = numRows;
    mNumColumns = numColumns;

    ReadAccess<IndexType> csrIA( ia, loc );
    ReadAccess<IndexType> csrJA( ja, loc );
    ReadAccess<OtherValueType> csrValues( values, loc );

    if( !OpenMPCSRUtils::validOffsets( csrIA.get(), numRows, numValues ) )
    {
        COMMON_THROWEXCEPTION( "invalid offset array" )
    }

    if( !OpenMPUtils::validIndexes( csrJA.get(), numValues, numColumns ) )
    {
        COMMON_THROWEXCEPTION( "invalid column indexes, #columns = " << numColumns )
    }

    WriteOnlyAccess<ValueType> data( mData, loc, mNumRows * mNumColumns );

    OpenMPDenseUtils::setCSRValues( data.get(), mNumRows, mNumColumns, csrIA.get(), csrJA.get(), csrValues.get() );

    mDiagonalProperty = checkDiagonalProperty();

    // dense storage does not care about diagonal property, is always okay
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::invert( const MatrixStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "invert( " << other << ") to a dense storage" )

    if( other.getFormat() == Format::DENSE )
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

    if( &other != this )
    {
        SCAI_LOG_INFO( logger, "invertDense: copy input matrix to this matrix" )
        assign( other );
    }

    int nRows = this->getNumRows();
    int nCols = this->getNumColumns();

    SCAI_ASSERT_EQUAL_ERROR( nRows, nCols )

    static LAMAKernel<BLASInterface::getinv<ValueType> > getinv;

    ContextPtr loc = getinv.getValidContext( this->getContextPtr() );

    WriteAccess<ValueType> denseValues( this->getData(), loc );

    SCAI_CONTEXT_ACCESS( loc );

    getinv[loc]( nRows, denseValues.get(), nCols );

    SCAI_LOG_INFO( logger, "invertDense: this = " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::matrixTimesVector(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y" << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    if( mNumRows == 0 )
    {
        return; // nothing to do
    }

    ContextPtr loc = mContext;

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    // not used here: LAMA_INTERFACE_FN_T( normalGEMV, loc, DenseUtils, Mult, ValueType )

    // using BLAS2 interface requires result and y to be aliased

    if( beta == scai::common::constants::ZERO )
    {
        SCAI_LOG_INFO( logger, "set result = 0 as y != result and beta = 0" )

        static LAMAKernel<UtilsInterface::setVal<ValueType> > setVal;

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), mNumRows, static_cast<ValueType>(0.0) );
    }
    else if( &result != &y )
    {
        SCAI_LOG_INFO( logger, "set result = y as y != result" )

        static LAMAKernel<BLASInterface::copy<ValueType> > copy;

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

        // by setting result = y we get the correct results

        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        copy[loc]( mNumRows, rY.get(), 1, wResult.get(), 1, NULL );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * A * x + beta * result

    if( mNumColumns == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if( beta == scai::common::constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if( beta == scai::common::constants::ONE )
        {
            // no scaling required
        }
        else
        {
            static LAMAKernel<BLASInterface::scal<ValueType> > scal;

            WriteAccess<ValueType> wResult( result, loc );

            SCAI_CONTEXT_ACCESS( loc )

            scal[loc]( mNumRows, beta, wResult.get(), 1, NULL );
        }
    }
    else
    {
        // mNumColums > 0, mnumRows > 0, so we avoid problems for gemv with m==0 or n==0

        static LAMAKernel<BLASInterface::gemv<ValueType> > gemv;

        ReadAccess<ValueType> denseValues( mData, loc );
        ReadAccess<ValueType> rX( x, loc );

        int lda = mNumColumns; // stride for denseValues between rows

        WriteAccess<ValueType> wResult( result, loc );

        SCAI_CONTEXT_ACCESS( loc )

        // gemv:  result = alpha * this * x + beta * result

        gemv[loc]( CblasRowMajor, CblasNoTrans, mNumRows, mNumColumns, alpha, denseValues.get(), lda, rX.get(),
                   1, beta, wResult.get(), 1, NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::vectorTimesMatrix(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y" << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )

    if( mNumRows == 0 )
    {
        return; // nothing to do
    }

    ContextPtr loc = mContext;

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    // not used here: LAMA_INTERFACE_FN_T( normalGEVM, loc, DenseUtils, Mult, ValueType )

    // using BLAS2 interface requires result and y to be aliased

    if( beta == scai::common::constants::ZERO )
    {
        SCAI_LOG_INFO( logger, "set result = 0 as y != result and beta = 0" )

        static LAMAKernel<UtilsInterface::setVal<ValueType> > setVal;

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wResult.get(), mNumColumns, static_cast<ValueType>(0.0) );
    }
    else if( &result != &y )
    {
        SCAI_LOG_INFO( logger, "set result = y as y != result" )

        static LAMAKernel<BLASInterface::copy<ValueType> > copy;

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

        // by setting result = y we get the correct results

        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        copy[loc]( mNumColumns, rY.get(), 1, wResult.get(), 1, NULL );
    }
    else
    {
        SCAI_LOG_INFO( logger, "alias of result and y, can use it" )
    }

    // now we have: result = alpha * A * x + beta * result

    if( mNumColumns == 0 )
    {
        SCAI_LOG_INFO( logger, "empty matrix, so compute result = " << beta << " * result " )

        if( beta == scai::common::constants::ZERO )
        {
            // nothing more to do, y is already 0
        }
        else if( beta == scai::common::constants::ONE )
        {
            // no scaling required
        }
        else
        {
            static LAMAKernel<BLASInterface::scal<ValueType> > scal;

            WriteAccess<ValueType> wResult( result, loc );

            SCAI_CONTEXT_ACCESS( loc )

            scal[loc]( mNumColumns, beta, wResult.get(), 1, NULL );
        }
    }
    else
    {
        // mNumColums > 0, mnumRows > 0, so we avoid problems for gevm with m==0 or n==0

        static LAMAKernel<BLASInterface::gemv<ValueType> > gemv;

        ReadAccess<ValueType> denseValues( mData, loc );
        ReadAccess<ValueType> rX( x, loc );

        int lda = mNumColumns; // stride for denseValues between rows

        WriteAccess<ValueType> wResult( result, loc );

        SCAI_CONTEXT_ACCESS( loc )

        // gemv:  result = alpha * this * x + beta * result

        gemv[loc]( CblasRowMajor, CblasTrans, mNumRows, mNumColumns, alpha, denseValues.get(), lda, rX.get(),
              1, beta, wResult.get(), 1, NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::print() const
{
    using std::cout;
    using std::endl;

    ReadAccess<ValueType> values( mData );

    cout << "DenseStorage " << mNumRows << " x " << mNumColumns << ", addr  = " << values.get() << endl;

    for( IndexType i = 0; i < mNumRows; i++ )
    {
        cout << "Row " << i << " :";

        for( IndexType j = 0; j < mNumColumns; j++ )
        {
            cout << " " << values[i * mNumColumns + j];
        }

        cout << endl;
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

    if( a.getFormat() == Format::DENSE )
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

    if( b.getFormat() == Format::DENSE )
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

    if( c.getFormat() == Format::DENSE )
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

    if( &a == this )
    {
        SCAI_LOG_INFO( logger, "temporary for A in A * B ( dense storages) needed" )
        tmpA = common::shared_ptr<DenseStorage<ValueType> >( new DenseStorage<ValueType>( a ) );
        ptrA = tmpA.get();
    }

    if( &b == this )
    {

        // avoid two temporaries

        if( &a == &b )
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

    ContextPtr context = mContext;

    if( beta == scai::common::constants::ZERO )
    {
        // do not care at all about C as it might be any dummy, or aliased to result

        static LAMAKernel<UtilsInterface::setVal<ValueType> > setVal;

        SCAI_LOG_INFO( logger, "init this result with 0, size = " << m * n )
        WriteOnlyAccess<ValueType> resAccess( getData(), context, m * n );
        SCAI_CONTEXT_ACCESS( context )
        setVal[context]( resAccess.get(), m * n, static_cast<ValueType>(0.0) );
    }
    else if( this != &c )
    {
        // force result = C

        SCAI_ASSERT_EQUAL_ERROR( m, c.getNumRows() )
        SCAI_ASSERT_EQUAL_ERROR( n, c.getNumColumns() )
        mNumRows = m;
        mNumColumns = n;
        static LAMAKernel<BLASInterface::copy<ValueType> > copy;
        ReadAccess<ValueType> cAccess( c.getData(), context );
        WriteOnlyAccess<ValueType> resAccess( getData(), context, m * n );
        SCAI_LOG_TRACE( logger, "Copying: res = c " )
        SCAI_CONTEXT_ACCESS( context )
        copy[context]( n * m, cAccess.get(), 1, resAccess.get(), 1, NULL );
    }
    else
    {
        SCAI_LOG_INFO( logger, "results is aliased with C as required for gemm, beta = " << beta )
    }

    // now C used for BLAS3 call of GEMM is this matrix

    ReadAccess<ValueType> aAccess( ptrA->getData(), context );
    ReadAccess<ValueType> bAccess( ptrB->getData(), context );
    WriteAccess<ValueType> resAccess( getData(), context );

    int lda = a.getNumColumns();
    int ldb = b.getNumColumns();
    int ldc = mNumColumns;

    // Now  C = alpha * A * B + beta * C, can use GEMM of BLAS3

    SCAI_LOG_TRACE( logger,
                    "GEMM( CblasRowMajor, CblasNoTrans, CblasNoTrans, m:" << m << ", n:" << n << ", k:" << k << ", alpha:" << alpha << ", A , lda:" << lda << " ,B ,ldb:"<< ldb << " ,beta:" << beta << " C ,ldc:" << ldc << " )" )

    if( lda != 0 && n != 0 && m != 0 )
    {
        static LAMAKernel<BLASInterface::gemm<ValueType> > gemm;

        SCAI_CONTEXT_ACCESS( context )

        gemm[context]( CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, aAccess.get(), lda, bAccess.get(), ldb, beta,
                       resAccess.get(), ldc, NULL );
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

    if( n == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<BLASInterface::asum<ValueType> > asum;

    ContextPtr loc = asum.getValidContext( this->getContextPtr() );

	ReadAccess<ValueType> data( mData, loc );

	SCAI_CONTEXT_ACCESS( loc );

    IndexType inc = 1;

	return asum[loc]( n, data.get(), inc, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorageView<ValueType>::l2Norm() const
{
    IndexType n = mNumRows * mNumColumns;

    if( n == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<BLASInterface::dot<ValueType> > dot;

    ContextPtr loc = dot.getValidContext( this->getContextPtr() );

	ReadAccess<ValueType> data( mData, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return ::sqrt(dot[loc]( n, data.get(), 1, data.get(), 1, NULL ));
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DenseStorageView<ValueType>::maxNorm() const
{
    IndexType n = mNumRows * mNumColumns;

    if( n == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<UtilsInterface::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = absMaxVal.getValidContext( this->getContextPtr() );

    ReadAccess<ValueType> read1( mData, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal[loc]( read1.get(), n );

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

    if( other.getValueType() == getValueType() && ( other.getFormat() == Format::DENSE ) )
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

    if( n == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<UtilsInterface::absMaxDiffVal<ValueType> > absMaxDiffVal;

    ContextPtr loc = absMaxDiffVal.getValidContext( this->getContextPtr() );

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

    // @ToDo: allow for arbitrary locations

    WriteOnlyAccess<ValueType> data( mData, mNumRows * mNumColumns );
    ReadAccess<OtherValueType> otherData( other.getData() );

    OpenMPDenseUtils::copyDenseValues( data.get(), mNumRows, mNumColumns, otherData.get() );

    SCAI_LOG_INFO( logger, *this << ": assigned dense storage " << other )

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DenseStorageView<ValueType>::assign( const _MatrixStorage& other )
{
    // Here more efficient solutions can be used instead of build/set CSR data

    if( &other == this )
    {
        SCAI_LOG_INFO( logger, "self assign, is skipped" )
        return;
    }

    if( other.getFormat() == Format::DENSE )
    {
        // more efficient solution for assigment of dense storage

        ScalarType arrayType = other.getValueType();

        switch( arrayType )
        {

#define LAMA_ASSIGN_DENSE_CALL( z, I, _ )                                                \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I:                                          \
{                                                                                        \
    const DenseStorageView<ARITHMETIC_HOST_TYPE_##I>* otherTyped =                       \
            dynamic_cast<const DenseStorageView<ARITHMETIC_HOST_TYPE_##I>*>( &other );   \
    SCAI_ASSERT_DEBUG( otherTyped, other << ": dynamic cast failed, should not happen" ) \
    assignDenseStorageImpl( *otherTyped );                                               \
    return;                                                                              \
}                                                                                        \

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_ASSIGN_DENSE_CALL, _ )

#undef LAMA_ASSIGN_DENSE_CALL

default            :

    SCAI_LOG_INFO( logger, "unsupported  typed assign" )

    // also take fallback call
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

ContextPtr loc = mContext;

WriteOnlyAccess<ValueType> data( mData, loc, numRows * numCols );

SCAI_LOG_DEBUG( logger, *this << " allocated, #values = " << mNumRows * mNumColumns << ", not initialized" )
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
                "get value (" << i << ", " << j << ")--> " << i*mNumColumns+j << "= " << data[i*mNumColumns+j] )
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

WriteOnlyAccess<ValueType> data( mData, numRows * numColumns );

OpenMPUtils::setVal( data.get(), numRows * numColumns, static_cast<ValueType>(0.0) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorage<ValueType>::DenseStorage(
    const LAMAArray<ValueType>& data,
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
DenseStorage<ValueType>::DenseStorage( const _MatrixStorage& other, const ContextPtr loc )

: DenseStorageView<ValueType>( mDataArray, 0, 0, false )

{
_MatrixStorage::setContext( loc );
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
DenseStorageView<ValueType>::clone() const
{
return new DenseStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DenseStorageView<ValueType>*
DenseStorageView<ValueType>::copy() const
{
return new DenseStorage<ValueType>( *this );
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_DENSE_STORAGE_INSTANTIATE(z, I, _)                                         \
    template<>                                                                          \
    const char* DenseStorage<ARITHMETIC_HOST_TYPE_##I>::typeName()                      \
    {                                                                                   \
        return "DenseStorage<" PRINT_STRING(ARITHMETIC_HOST_TYPE_##I) ">";                \
    }                                                                                   \
                                                                                        \
    template<>                                                                          \
    const char* DenseStorageView<ARITHMETIC_HOST_TYPE_##I>::typeName()                  \
    {                                                                                   \
    	return "DenseStorage<" PRINT_STRING(ARITHMETIC_HOST_TYPE_##I) ">";                \
    }                                                                                   \
                                                                                        \
    template class COMMON_DLL_IMPORTEXPORT DenseStorage<ARITHMETIC_HOST_TYPE_##I> ;     \
    template class COMMON_DLL_IMPORTEXPORT DenseStorageView<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_DENSE_STORAGE_INSTANTIATE, _ )

#undef LAMA_DENSE_STORAGE_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
