/**
 * @file COOStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
#include <lama/storage/COOStorage.hpp>

// others
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>
#include <lama/ContextFactory.hpp>
#include <lama/ContextAccess.hpp>
#include <lama/LAMAInterface.hpp>

#include <lama/LAMAArrayUtils.hpp>

#include <lama/openmp/OpenMPUtils.hpp>
#include <lama/openmp/OpenMPCOOUtils.hpp>
#include <lama/openmp/OpenMPCSRUtils.hpp>

#include <lama/tracing.hpp>

using std::auto_ptr;
using boost::shared_ptr;

namespace lama
{

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, COOStorage<ValueType>::logger, "MatrixStorage.COOStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage(const IndexType numRows, const IndexType numColumns) :

    CRTPMatrixStorage<COOStorage<ValueType>, ValueType> ( numRows, numColumns ),
    mNumValues(0)
{
    LAMA_LOG_DEBUG( logger, "COOStorage for matrix " << mNumRows << " x " << mNumColumns
                    << ", no non-zero elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::COOStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<ValueType>& values )

    : CRTPMatrixStorage<COOStorage<ValueType>,ValueType>( )
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
    LAMA_LOG_DEBUG( logger, "COOStorage, matrix is 0 x 0." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat COOStorage<ValueType>::getFormat() const
{
    return COO;
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

        HostReadAccess<IndexType> ia( mIA );
        HostReadAccess<IndexType> ja( mJA );

        // The diagonal property is given if the first numDiags entries
        // are the diagonal elements

        #pragma omp parallel for schedule(LAMA_OMP_SCHEDULE)
        for ( IndexType i = 0; i < mNumRows; ++i )
        {
            if ( !diagonalProperty )
            {
                continue;
            }
    
            if ( ia[i] != i || ja[i] != i )
            {
                diagonalProperty = false;
            }
        }
    }

    LAMA_LOG_INFO( logger, *this << ": checkDiagonalProperty -> " << diagonalProperty )
   
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
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mIA.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mJA.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )

    // check row indexes in IA and column indexes in JA

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT( validIndexes, loc, Utils, Indexes )

        ReadAccess<IndexType> rJA( mJA, loc );
        ReadAccess<IndexType> rIA( mIA, loc );

        LAMA_CONTEXT_ACCESS( loc )

        LAMA_ASSERT_ERROR( validIndexes ( rIA.get(), mNumValues, mNumRows ),
                           *this << " @ " << msg << ": illegel indexes in IA" )

        LAMA_ASSERT_ERROR( validIndexes ( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setIdentity( const IndexType size )
{
    LAMA_LOG_INFO( logger, "set identity values for " << size )

    mNumRows = size;
    mNumColumns = size;
    mNumValues = mNumRows;

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )
    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

    WriteOnlyAccess<IndexType> ia( mIA, loc, mNumValues );
    WriteOnlyAccess<IndexType> ja( mJA, loc, mNumValues );
    WriteOnlyAccess<ValueType> values( mValues, loc, mNumValues );

    ValueType one = static_cast<ValueType>( 1.0 );

    LAMA_CONTEXT_ACCESS( loc )

    setOrder( ia.get(), mNumValues );
    setOrder( ja.get(), mNumValues );

    setVal( values.get(), mNumValues, one );

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

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    LAMA_INTERFACE_FN( getCSRSizes, loc, COOUtils, Counting )
    LAMA_INTERFACE_FN( sizes2offsets, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_TT( getCSRValues, loc, COOUtils, Conversions, ValueType, OtherValueType )

    WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );
    ReadAccess<IndexType> cooIA( mIA, loc );

    getCSRSizes( csrIA.get(), mNumRows, mNumValues, cooIA.get() );

    if ( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets( csrIA.get(), mNumRows );

    LAMA_ASSERT_EQUAL_DEBUG( mNumValues, numValues )

    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );

    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<OtherValueType> csrValues( *values, loc, numValues );

    getCSRValues( csrJA.get(), csrValues.get(), csrIA.get(), mNumRows, mNumValues, cooIA.get(),
                  cooJA.get(), cooValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::setCOOData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const _LAMAArray& values )
{
    // check the sizes of the arrays

    LAMA_ASSERT_EQUAL_ERROR( numValues, ia.size() )
    LAMA_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    LAMA_ASSERT_EQUAL_ERROR( numValues, values.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValues = numValues;

    ContextPtr loc = getContextPtr();

    LAMAArrayUtils::assignImpl( mIA, ia, loc );
    LAMAArrayUtils::assignImpl( mJA, ja, loc );

    LAMAArrayUtils::assign( mValues, values, loc );  // supports type conversion

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode

#ifdef LAMA_ASSERT_LEVEL_DEBUG
    check( "COOStorage.setCOOData" );
#endif

    mDiagonalProperty = checkDiagonalProperty();

    // Note: no support for row indexes in COO format

    LAMA_LOG_INFO( logger, *this << ": set COO by arrays ia, ja, values" )
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
    ContextPtr loc = ContextFactory::getContext( Context::Host );

    ReadAccess<IndexType> csrJA( ja, loc );
    ReadAccess<OtherValueType> csrValues( values, loc );

    mNumRows = numRows;
    mNumColumns = numColumns;

    // check if input csr data has the diagonal property and inherit it

    int numDiagonals = std::min( numRows, numColumns );

    {
        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<IndexType> csrJA( ja, loc );
        mDiagonalProperty = OpenMPCSRUtils::hasDiagonalProperty( numDiagonals, csrIA.get(), csrJA.get() );
    }

    if ( !mDiagonalProperty )
    {
        numDiagonals = 0; // do not store diagonal data at the beginning in COO data
    }

    mNumValues = numValues;

    LAMA_LOG_DEBUG( logger,
                    "input csr data with " << mNumValues << "entries,  has diagonal property = " << mDiagonalProperty )

    {
        LAMA_INTERFACE_FN( offsets2ia, loc, COOUtils, Counting );

        ReadAccess<IndexType> csrIA( ia, loc );
        WriteOnlyAccess<IndexType> cooIA( mIA, loc, mNumValues );
        offsets2ia( cooIA.get(), mNumValues, csrIA.get(), mNumRows, numDiagonals );
    }

    {
        LAMA_INTERFACE_FN_TT( setCSRData, loc, COOUtils, Conversions, IndexType, IndexType );

        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<IndexType> csrJA( ja, loc );
        WriteOnlyAccess<IndexType> cooJA( mJA, loc, mNumValues );

        setCSRData( cooJA.get(), csrJA.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }

    {
        LAMA_INTERFACE_FN_TT( setCSRData, loc, COOUtils, Conversions, ValueType, OtherValueType );

        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<OtherValueType> csrValues( values, loc );
        WriteOnlyAccess<ValueType> cooValues( mValues, loc, mNumValues );

        setCSRData( cooValues.get(), csrValues.get(), numValues, csrIA.get(), mNumRows, numDiagonals );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>::~COOStorage()
{
    LAMA_LOG_DEBUG( logger, "~COOStorage for matrix " << mNumRows << " x " << mNumColumns )
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
    LAMA_LOG_INFO( logger, "allocate COO sparse matrix of size " << numRows << " x " << numColumns )

    clear();   // all variables are set for a zero-sized matrix

    mNumRows = numRows;
    mNumColumns = numColumns;

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "COO(rows=" << mNumRows << ",cols=" << mNumColumns << ")";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType COOStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    // only supported on Host at this time

    const HostReadAccess<IndexType> ia( mIA );
    const HostReadAccess<IndexType> ja( mJA );
    const HostReadAccess<ValueType> values( mValues );

    LAMA_LOG_DEBUG( logger, "get value (" << i << ", " << j << ") from " << *this )

    for ( IndexType kk = 0; kk < mNumValues; ++kk )
    {
        if ( ia[kk] == i && ja[kk] == j )
        {
            return values[kk];
        }
    }

    return 0.0;
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

    HostWriteAccess<ValueType> wValues( mValues );
    HostReadAccess<IndexType> rJa( mJA );
    HostReadAccess<IndexType> rIa( mIA );
    ValueType value = scalar.getValue<ValueType>();

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wValues[i] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    HostWriteAccess<ValueType> wValues( mValues );
    ValueType value = scalar.getValue<ValueType>();

    for ( IndexType i = 0; i < mNumValues; ++i )
    {
        wValues[i] *= value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::scaleImpl( const LAMAArray<OtherType>& values )
{
    HostReadAccess<OtherType> rValues( values );
    HostWriteAccess<ValueType> wValues( mValues );
    HostReadAccess<IndexType> rIa( mIA );

    for ( IndexType i = 0; i < mNumValues; ++i )
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
    LAMA_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    HostWriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const HostReadAccess<IndexType> ia( mIA );
    const HostReadAccess<IndexType> ja( mJA );
    const HostReadAccess<ValueType> values( mValues );

    for ( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = 0.0;
    }

    for ( IndexType kk = 0; kk < mNumValues; ++kk )
    {
        if ( ia[kk] != i )
        {
            continue;
        }

        wRow[ja[kk]] = static_cast<OtherType>( values[kk] );
    }
}

/* --------------------------------------------------------------------------- */

// Note: template instantation of this method for OtherType=[double,float] is
//       done implicitly by getDiagonal method of CRTPMatrixStorage

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, OtherType, ValueType )

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    // diagonal elements are the first entries of mValues

    set( wDiagonal.get(), rValues.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

// Note: template instantation of this method for OtherType=[double,float] is
//       done implicitly by setDiagonal method of CRTPMatrixStorage

template<typename ValueType>
template<typename OtherType>
void COOStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherType>& diagonal )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, ValueType, OtherType )

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    // diagonal elements are the first entries of mValues

    set( wValues.get(), rDiagonal.get(), numDiagonalElements );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType COOStorage<ValueType>::maxNorm() const
{
    LAMA_LOG_INFO( logger, *this << ": maxNorm()" )

    const IndexType n = mNumValues;

    if ( n == 0 )
    {
        return 0.0f;
    }

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, Utils, Reductions, ValueType )

    ReadAccess<ValueType> cooValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal( cooValues.get(), n );

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

    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const

{
    LAMA_REGION( "Storage.COO.timesVector" )

    LAMA_LOG_DEBUG( logger,
                    "Computing z = alpha * A * x + beta * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    // Method on CUDA is not safe due to atomic

    ContextPtr loc = getContextPtr();

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_DEFAULT_T( normalGEMV, loc, COOUtils, Mult, ValueType )

    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        LAMA_CONTEXT_ACCESS( loc )
        normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, cooIA.get(), cooJA.get(),
                    cooValues.get(), mNumValues, NULL );
    }
    else
    {
        // make also sure that result will have the correct size

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc )
        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, cooIA.get(), cooJA.get(), cooValues.get(),
                    mNumValues, NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
auto_ptr<SyncToken> COOStorage<ValueType>::matrixTimesVectorAsyncToDo(
    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const
{
    LAMA_LOG_DEBUG( logger,
                    "Computing z = alpha * A * x + beta * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    // not yet available on other devices, so we take Host

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVectorAsync on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEMV, loc, COOUtils, Mult, ValueType )

    auto_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > cooIA( new ReadAccess<IndexType>( mIA, loc ) );
    shared_ptr<ReadAccess<IndexType> > cooJA( new ReadAccess<IndexType>( mJA, loc ) );
    shared_ptr<ReadAccess<ValueType> > cooValues( new ReadAccess<ValueType>( mValues, loc ) );
    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        // we assume that normalGEMV can deal with the alias of result, y

        LAMA_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, cooIA->get(), cooJA->get(),
                    cooValues->get(), mNumValues, syncToken.get() );

        syncToken->pushAccess( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumRows ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        LAMA_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, cooIA->get(), cooJA->get(),
                    cooValues->get(), mNumValues, syncToken.get() );

        syncToken->pushAccess( shared_ptr<BaseAccess>( wResult ) );
        syncToken->pushAccess( shared_ptr<BaseAccess>( rY ) );
    }

    syncToken->pushAccess( cooIA );
    syncToken->pushAccess( cooJA );
    syncToken->pushAccess( cooValues );
    syncToken->pushAccess( rX );

    return syncToken;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void COOStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.COO.jacobiIterate" )

    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    LAMA_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( solution == oldSolution )
    {
        LAMA_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square

    ContextPtr loc = getContextPtr();

    // loc = ContextFactory::getContext( Context::Host );  // does not run on other devices

    LAMA_INTERFACE_FN_DEFAULT_T( jacobi, loc, COOUtils, Solver, ValueType )

    WriteAccess<ValueType> wSolution( solution, loc );
    ReadAccess<IndexType> cooIA( mIA, loc );
    ReadAccess<IndexType> cooJA( mJA, loc );
    ReadAccess<ValueType> cooValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    // Due to diagonal property there is no advantage by taking row indexes

    LAMA_CONTEXT_ACCESS( loc )

    jacobi( wSolution.get(), mNumValues, cooIA.get(), cooJA.get(), cooValues.get(), rOldSolution.get(), rRhs.get(),
            omega, mNumRows, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::create() const
{
    return new COOStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
COOStorage<ValueType>* COOStorage<ValueType>::copy() const
{
    return new COOStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<>
const char* COOStorage<float>::typeName()
{
    return "COOStorage<float>";
}

template class LAMA_DLL_IMPORTEXPORT COOStorage<float> ;

template<>
const char* COOStorage<double>::typeName()
{
    return "COOStorage<double>";
}

template class LAMA_DLL_IMPORTEXPORT COOStorage<double> ;

} // namespace lama
