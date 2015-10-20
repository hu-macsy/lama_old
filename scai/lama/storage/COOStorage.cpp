/**
 * @file COOStorage.cpp
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
 * @brief Instantitions for template class COOStorage.
 * @author Lauretta Schubert
 * @date 25.05.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/storage/COOStorage.hpp>

// local library
#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/LAMAArrayUtils.hpp>
#include <scai/lama/kernel_registry.hpp>

#include <scai/lama/openmp/OpenMPUtils.hpp>
#include <scai/lama/openmp/OpenMPCOOUtils.hpp>
#include <scai/lama/openmp/OpenMPCSRUtils.hpp>

// internal scai libraries
#include <scai/hmemo.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>

// boost
#include <boost/preprocessor.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::unique_ptr;
using common::shared_ptr;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, COOStorage<ValueType>::logger, "MatrixStorage.COOStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( const IndexType numRows, const IndexType numColumns )
    :

    CRTPMatrixStorage<COOStorage<ValueType>,ValueType>( numRows, numColumns ), mNumValues( 0 )
{
    SCAI_LOG_DEBUG( logger, "COOStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<ValueType>& values )

    : CRTPMatrixStorage<COOStorage<ValueType>,ValueType>()
{
    // all array must have the same size

    IndexType numValues = ia.size();

    setCOOData( numRows, numColumns, numValues, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage( const COOStorage<ValueType>& other )

    : CRTPMatrixStorage<COOStorage<ValueType>,ValueType>( 0, 0 )
{
    // ToDo: copy of same storage format should be more efficient

    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage()
    : CRTPMatrixStorage<COOStorage<ValueType>,ValueType>( 0, 0 ), mNumValues( 0 )
{
    SCAI_LOG_DEBUG( logger, "COOStorage, matrix is 0 x 0." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat COOStorage<ValueType>::getFormat() const
{
    return Format::COO;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool COOStorage<ValueType>::checkDiagonalProperty() const
{
    bool diagonalProperty = true;

    if( mNumRows != mNumColumns )
    {
        diagonalProperty = false;
    }
    else if( mNumRows == 0 )
    {
        // zero sized matrix

        diagonalProperty = true;
    }
    else if( mIA.size() == 0 )
    {
        diagonalProperty = false;
    }
    else
    {
        diagonalProperty = true; // intialization for reduction

        ContextPtr contextPtr = Context::getHostPtr();

        ReadAccess<IndexType> ia( mIA, contextPtr );
        ReadAccess<IndexType> ja( mJA, contextPtr );

        // The diagonal property is given if the first numDiags entries
        // are the diagonal elements

        #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

        for( IndexType i = 0; i < mNumRows; ++i )
        {
            if( !diagonalProperty )
            {
                continue;
            }

            if( ia[i] != i || ja[i] != i )
            {
                diagonalProperty = false;
            }
        }
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
        scai::kregistry::KernelTraitContextFunction<UtilsInterface::validIndexes> validIndexes;

        ContextPtr  loc = getValidContext( getContextPtr(), validIndexes );  // find location where routine is available
        ContextType ctx = loc->getType();

        ReadAccess<IndexType> rJA( mJA, loc );
        ReadAccess<IndexType> rIA( mIA, loc );

        SCAI_CONTEXT_ACCESS( loc )

        SCAI_ASSERT_ERROR( validIndexes[ ctx ]( rIA.get(), mNumValues, mNumRows ),
                           *this << " @ " << msg << ": illegel indexes in IA" )

        SCAI_ASSERT_ERROR( validIndexes[ ctx ]( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
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

    scai::kregistry::KernelTraitContextFunction<UtilsInterface::setOrder<IndexType> > setOrder;
    scai::kregistry::KernelTraitContextFunction<UtilsInterface::setVal<ValueType> > setVal;

    ContextPtr loc = getValidContext( getContextPtr(), setOrder );

    common::ContextType ctx = loc->getType();

    WriteOnlyAccess<IndexType> ia( mIA, loc, mNumValues );
    WriteOnlyAccess<IndexType> ja( mJA, loc, mNumValues );
    WriteOnlyAccess<ValueType> values( mValues, loc, mNumValues );

    SCAI_CONTEXT_ACCESS( loc )

    setOrder[ ctx ]( ia.get(), mNumValues );
    setOrder[ ctx ]( ja.get(), mNumValues );

    setVal[ ctx ]( values.get(), mNumValues, static_cast<ValueType>(1.0) );

    mDiagonalProperty = true;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::buildCSR(
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr /* loc */) const
{
    // multiple routines from interface needed, so do it on Host to be safe

    ContextPtr loc = Context::getContextPtr( context::Host );

    LAMA_INTERFACE_FN( getCSRSizes, loc, COOUtils, Counting )
    LAMA_INTERFACE_FN( sizes2offsets, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_TT( getCSRValues, loc, COOUtils, Conversions, ValueType, OtherValueType )

    WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );
    ReadAccess<IndexType> cooIA( mIA, loc );

    getCSRSizes( csrIA.get(), mNumRows, mNumValues, cooIA.get() );

    if( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets( csrIA.get(), mNumRows );

    SCAI_ASSERT_EQUAL_DEBUG( mNumValues, numValues )

    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );

    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<OtherValueType> csrValues( *values, loc, numValues );

    getCSRValues( csrJA.get(), csrValues.get(), csrIA.get(), mNumRows, mNumValues, cooIA.get(), cooJA.get(),
                  cooValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setCOOData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const ContextArray& values )
{
    // check the sizes of the arrays

    SCAI_ASSERT_EQUAL_ERROR( numValues, ia.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, values.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValues = numValues;

    ContextPtr loc = getContextPtr();

    LAMAArrayUtils::assignImpl( mIA, ia, loc );
    LAMAArrayUtils::assignImpl( mJA, ja, loc );

    LAMAArrayUtils::assign( mValues, values, loc ); // supports type conversion

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode

#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "COOStorage.setCOOData" );
#endif

    mDiagonalProperty = checkDiagonalProperty();

    // Note: no support for row indexes in COO format

    SCAI_LOG_INFO( logger, *this << ": set COO by arrays ia, ja, values" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void COOStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    const ContextPtr )
{
    SCAI_LOG_DEBUG( logger, "set CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues )

    SCAI_ASSERT_EQUAL_DEBUG( numRows + 1, ia.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, values.size() )

    ContextPtr loc = getContextPtr();

    ReadAccess<IndexType> csrJA( ja, loc );
    ReadAccess<OtherValueType> csrValues( values, loc );

    mNumRows = numRows;
    mNumColumns = numColumns;

    // check if input csr data has the diagonal property and inherit it

    int numDiagonals = std::min( numRows, numColumns );

    {
        SCAI_LOG_DEBUG( logger,
                        "check CSR data " << numRows << " x " << numColumns << ", nnz = " << numValues << " for diagonal property, #diagonals = " << numDiagonals )

        LAMA_INTERFACE_FN( hasDiagonalProperty, loc, CSRUtils, Offsets );

        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<IndexType> csrJA( ja, loc );

        SCAI_CONTEXT_ACCESS( loc )

        mDiagonalProperty = hasDiagonalProperty( numDiagonals, csrIA.get(), csrJA.get() );
    }

    if( !mDiagonalProperty )
    {
        numDiagonals = 0; // do not store diagonal data at the beginning in COO data
    }

    mNumValues = numValues;

    SCAI_LOG_DEBUG( logger,
                    "input csr data with " << mNumValues << "entries,  has diagonal property = " << mDiagonalProperty )

    {
        LAMA_INTERFACE_FN( offsets2ia, loc, COOUtils, Counting );

        ReadAccess<IndexType> csrIA( ia, loc );
        WriteOnlyAccess<IndexType> cooIA( mIA, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        offsets2ia( cooIA.get(), mNumValues, csrIA.get(), mNumRows, numDiagonals );
    }

    {
        LAMA_INTERFACE_FN_TT( setCSRData, loc, COOUtils, Conversions, IndexType, IndexType );

        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<IndexType> csrJA( ja, loc );
        WriteOnlyAccess<IndexType> cooJA( mJA, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        setCSRData( cooJA.get(), csrJA.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }

    {
        LAMA_INTERFACE_FN_TT( setCSRData, loc, COOUtils, Conversions, ValueType, OtherValueType );

        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<OtherValueType> csrValues( values, loc );
        WriteOnlyAccess<ValueType> cooValues( mValues, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        setCSRData( cooValues.get(), csrValues.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }
}

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
    // only supported on Host at this time

    ContextPtr loc = Context::getHostPtr();

    const ReadAccess<IndexType> ia( mIA, loc );
    const ReadAccess<IndexType> ja( mJA, loc );
    const ReadAccess<ValueType> values( mValues, loc );

    SCAI_LOG_DEBUG( logger, "get value (" << i << ", " << j << ") from " << *this )

    for( IndexType kk = 0; kk < mNumValues; ++kk )
    {
        if( ia[kk] == i && ja[kk] == j )
        {
            return values[kk];
        }
    }

    return static_cast<ValueType>(0.0);
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
const LAMAArray<IndexType>& COOStorage<ValueType>::getIA() const
{
    return mIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& COOStorage<ValueType>::getJA() const
{
    return mJA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<ValueType>& COOStorage<ValueType>::getValues() const
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
void COOStorage<ValueType>::setDiagonalImpl( const Scalar scalar )
{
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ContextPtr loc = Context::getHostPtr();

    WriteAccess<ValueType> wValues( mValues, loc );
    ReadAccess<IndexType> rJa( mJA, loc );
    ReadAccess<IndexType> rIa( mIA, loc );

    ValueType value = scalar.getValue<ValueType>();

    for( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wValues[i] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    WriteAccess<ValueType> wValues( mValues );
    ValueType value = scalar.getValue<ValueType>();

    for( IndexType i = 0; i < mNumValues; ++i )
    {
        wValues[i] *= value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::scaleImpl( const LAMAArray<OtherType>& values )
{
    ContextPtr loc = Context::getHostPtr();

    ReadAccess<OtherType> rValues( values, loc );
    WriteAccess<ValueType> wValues( mValues, loc );  // update
    ReadAccess<IndexType> rIa( mIA, loc );

    // Only host implementation available

    for( IndexType i = 0; i < mNumValues; ++i )
    {
        wValues[i] *= static_cast<ValueType>( rValues[rIa[i]] );
    }
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
    std::swap( mNumValues, other.mNumValues );
    mIA.swap( other.mIA );
    mJA.swap( other.mJA );
    mValues.swap( other.mValues );

    MatrixStorage<ValueType>::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    ContextPtr hostContext = Context::getHostPtr();

    WriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const ReadAccess<IndexType> ia( mIA, hostContext );
    const ReadAccess<IndexType> ja( mJA, hostContext );
    const ReadAccess<ValueType> values( mValues, hostContext );

    // ToDo: OpenMP parallelization, interface

    for( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = static_cast<OtherType>(0.0);
    }

    for( IndexType kk = 0; kk < mNumValues; ++kk )
    {
        if( ia[kk] != i )
        {
            continue;
        }

        wRow[ja[kk]] = static_cast<OtherType>( values[kk] );
    }
}

/* --------------------------------------------------------------------------- */

// Note: template instantation of this method for OtherType is
//       done implicitly by getDiagonal method of CRTPMatrixStorage
template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static kregistry::KernelTraitContextFunction<UtilsInterface::set<OtherType, ValueType> > set;

    ContextPtr loc = getValidContext( this->getContextPtr(), set );

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // diagonal elements are the first entries of mValues

    set[ loc->getType() ]( wDiagonal.get(), rValues.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

// Note: template instantation of this method for OtherType is
//       done implicitly by setDiagonal method of CRTPMatrixStorage
template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherType>& diagonal )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static kregistry::KernelTraitContextFunction<UtilsInterface::set<ValueType, OtherType> > set;

    ContextPtr loc = getValidContext( this->getContextPtr(), set );

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // diagonal elements are the first entries of mValues

    set[ loc->getType() ]( wValues.get(), rDiagonal.get(), numDiagonalElements );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType COOStorage<ValueType>::l1Norm() const
{
	SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    const IndexType n = mNumValues;

	ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( asum, loc, BLAS, BLAS1, ValueType )

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return asum( n, data.get(), static_cast<IndexType>(1.0), NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::l2Norm() const
{
	SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    const IndexType n = mNumValues;

	ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( dot, loc, BLAS, BLAS1, ValueType )

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return ::sqrt(dot( n, data.get(), 1, data.get(), 1, NULL ));
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    const IndexType n = mNumValues;

    if( n == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static kregistry::KernelTraitContextFunction<UtilsInterface::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = getValidContext( this->getContextPtr(), absMaxVal );

    ReadAccess<ValueType> cooValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal[ loc->getType() ]( cooValues.get(), n );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t COOStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType) * mIA.size();
    memoryUsage += sizeof(IndexType) * mJA.size();
    memoryUsage += sizeof(ValueType) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::matrixTimesVector(

    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const

{
    SCAI_REGION( "Storage.COO.timesVector" )

    SCAI_LOG_DEBUG( logger,
                    "Computing z = alpha * A * x + beta * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    // Method on CUDA is not safe due to atomic

    ContextPtr loc = getContextPtr();

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_DEFAULT_T( normalGEMV, loc, COOUtils, Mult, ValueType )

    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )
        normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumValues, cooIA.get(), cooJA.get(),
                    cooValues.get(), NULL );
    }
    else
    {
        // make also sure that result will have the correct size

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )
        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumValues, cooIA.get(), cooJA.get(),
                    cooValues.get(), NULL );
    }
}
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::vectorTimesMatrix(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.COO.VectorTimesMatrix" )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    ContextPtr loc = getContextPtr();

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEVM, loc, COOUtils, Mult, ValueType )

    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEVM can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )
        normalGEVM( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumColumns, mNumValues, cooIA.get(),
                    cooJA.get(), cooValues.get(), NULL );
    }
    else
    {
        // make also sure that result will have the correct size

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )
        normalGEVM( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumColumns, mNumValues, cooIA.get(), cooJA.get(),
                    cooValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::matrixTimesVectorAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_DEBUG( logger,
                    "Computing z = alpha * A * x + beta * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    // not yet available on other devices, so we take Host

    ContextPtr loc = Context::getContextPtr( context::Host );

    if ( loc->getType() == context::Host )
    {
        // execution as separate thread

        void (COOStorage::*pf)(
            LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>& ) const

            = &COOStorage<ValueType>::matrixTimesVector;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        SCAI_LOG_INFO( logger, *this << ": matrixTimesVectorAsync on Host by own thread" )

        return new TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVectorAsync on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEMV, loc, COOUtils, Mult, ValueType )

    unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > cooIA( new ReadAccess<IndexType>( mIA, loc ) );
    shared_ptr<ReadAccess<IndexType> > cooJA( new ReadAccess<IndexType>( mJA, loc ) );
    shared_ptr<ReadAccess<ValueType> > cooValues( new ReadAccess<ValueType>( mValues, loc ) );
    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumValues, cooIA->get(),
                    cooJA->get(), cooValues->get(), syncToken.get() );

        syncToken->pushToken( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumRows ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumValues, cooIA->get(), cooJA->get(),
                    cooValues->get(), syncToken.get() );

        syncToken->pushToken( wResult );
        syncToken->pushToken( rY );
    }

    syncToken->pushToken( cooIA );
    syncToken->pushToken( cooJA );
    syncToken->pushToken( cooValues );
    syncToken->pushToken( rX );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOStorage<ValueType>::vectorTimesMatrixAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.COO.vectorTimesMatrixAsync" )

    ContextPtr loc = getContextPtr();

    // Note: checks will be done by asynchronous task in any case
    //       and exception in tasks are handled correctly

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on " << *loc )

    if( loc->getType() == context::Host )
    {
        // execution as separate thread

        void (COOStorage::*pf)(
            LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>& ) const

            = &COOStorage<ValueType>::vectorTimesMatrix;

        SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on Host by own thread" )

        using scai::common::bind;
        using scai::common::ref;

        return new TaskSyncToken( bind( pf, this, ref( result ), alpha, ref( x ), beta, ref( y ) ) );
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    LAMA_INTERFACE_FN_T( normalGEVM, loc, COOUtils, Mult, ValueType )

    unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > cooIA( new ReadAccess<IndexType>( mIA, loc ) );
    shared_ptr<ReadAccess<IndexType> > cooJA( new ReadAccess<IndexType>( mJA, loc ) );
    shared_ptr<ReadAccess<ValueType> > cooValues( new ReadAccess<ValueType>( mValues, loc ) );
    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumValues, cooIA->get(),
                    cooJA->get(), cooValues->get(), syncToken.get() );

        syncToken->pushToken( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumColumns ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumValues, cooIA->get(), cooJA->get(),
                    cooValues->get(), syncToken.get() );

        syncToken->pushToken( wResult );
        syncToken->pushToken( rY );
    }

    syncToken->pushToken( cooIA );
    syncToken->pushToken( cooJA );
    syncToken->pushToken( cooValues );
    syncToken->pushToken( rX );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.COO.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square

    ContextPtr loc = getContextPtr();

    // loc = Context::getContextPtr( context::Host );  // does not run on other devices

    LAMA_INTERFACE_FN_DEFAULT_T( jacobi, loc, COOUtils, Solver, ValueType )

    WriteAccess<ValueType> wSolution( solution, loc );
    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    // Due to diagonal property there is no advantage by taking row indexes

    SCAI_CONTEXT_ACCESS( loc )

    jacobi( wSolution.get(), mNumValues, cooIA.get(), cooJA.get(), cooValues.get(), rOldSolution.get(), rRhs.get(),
            omega, mNumRows, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::clone() const
{
    return new COOStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::copy() const
{
    return new COOStorage<ValueType>( *this );
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

#define LAMA_COO_STORAGE_INSTANTIATE(z, I, _)                                     \
    template<>                                                                    \
    const char* COOStorage<ARITHMETIC_HOST_TYPE_##I>::typeName()                  \
    {                                                                             \
        return "COOStorage<" PRINT_STRING(ARITHMETIC_HOST_TYPE_##I) ">";      \
    }                                                                             \
                                                                                  \
    template class COMMON_DLL_IMPORTEXPORT COOStorage<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_COO_STORAGE_INSTANTIATE, _ )

#undef LAMA_COO_STORAGE_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
