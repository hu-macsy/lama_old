/**
 * @file DIAStorage.cpp
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
 * @brief Instantiation for template class DIAStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/storage/DIAStorage.hpp>

// local library
#include <scai/lama/openmp/OpenMPCSRUtils.hpp>
#include <scai/lama/openmp/OpenMPDIAUtils.hpp>

#include <scai/lama/LAMAInterface.hpp>

// internal scai libraries
#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Constants.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::scoped_array;
// Allow for shared_ptr<ValueType> instead of common::shared_ptr<ValueType>
using common::shared_ptr;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DIAStorage<ValueType>::logger, "MatrixStorage.DIAStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const IndexType numRows, const IndexType numColumns )

    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( numRows, numColumns ), mNumDiagonals( 0 )
{
    SCAI_LOG_DEBUG( logger, "DIAStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage()
    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( 0, 0 ), mNumDiagonals( 0 )
{
    SCAI_LOG_DEBUG( logger, "DIAStorage, matrix is 0 x 0." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const LAMAArray<IndexType>& offsets,
    const LAMAArray<ValueType>& values )

    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( numRows, numColumns ), mNumDiagonals(
          numDiagonals ), mOffset( offsets ), mValues( values )
{
    // set diagonal property inherited as given

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const DIAStorage<ValueType>& other )

    : CRTPMatrixStorage<DIAStorage<ValueType>,ValueType>( 0, 0 )
{
    // @todo: copy of same storage format should be implemented more efficiently

    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::print() const
{
    using std::cout;
    using std::endl;

    cout << "DIAStorage " << mNumRows << " x " << mNumColumns << ", #diags = " << mNumDiagonals << endl;

    ReadAccess<IndexType> offset( mOffset );
    ReadAccess<ValueType> values( mValues );

    cout << "Diagonal offsets:";

    for( IndexType d = 0; d < mNumDiagonals; d++ )
    {
        cout << " " << offset[d];
    }

    cout << endl;

    for( IndexType i = 0; i < mNumRows; i++ )
    {
        cout << "Row " << i << " :";

        for( IndexType ii = 0; ii < mNumDiagonals; ++ii )
        {
            const IndexType j = i + offset[ii];

            if( j < 0 )
            {
                continue;
            }

            if( j >= mNumColumns )
            {
                break;
            }

            cout << " " << j << ":" << values[i + ii * mNumRows];
        }

        cout << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;
    mNumDiagonals = 0;

    mOffset.clear();
    mValues.clear();

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat DIAStorage<ValueType>::getFormat() const
{
    return Format::DIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setDiagonalImpl( const Scalar scalar )
{
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": has not diagonal property, cannot set diagonal" )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ContextPtr loc = getContextPtr(); // take context of this storage to set

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

    {
        // not all values might be changed, so use WriteAccess instead of WriteOnlyAccess

        WriteAccess<ValueType> wValues( mValues, loc );
        ReadAccess<IndexType> rOffset( mOffset, loc );

        ValueType value = scalar.getValue<ValueType>();

        SCAI_CONTEXT_ACCESS( loc )

        setVal( wValues.get(), numDiagonalElements, value );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    WriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const ReadAccess<IndexType> offset( mOffset );
    const ReadAccess<ValueType> values( mValues );

    for( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = static_cast<OtherType>(0.0);
    }

    for( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if( j < 0 || j >= mNumColumns )
        {
            continue;
        }

        wRow[j] = static_cast<OtherType>( values[diaindex( i, d, mNumRows, mNumDiagonals )] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const
{
    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, OtherType, ValueType )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // Diagonal is first column

    set( wDiagonal.get(), rValues.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherType>& diagonal )
{
    ContextPtr loc = getContextPtr();

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    numDiagonalElements = std::min( numDiagonalElements, diagonal.size() );

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, ValueType, OtherType )

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    set( wValues.get(), rDiagonal.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    WriteAccess<ValueType> wValues( mValues );

    ValueType value = scalar.getValue<ValueType>();

    for( IndexType i = 0; i < mValues.size(); i++ )
    {
        wValues[i] *= value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::scaleImpl( const LAMAArray<OtherType>& diagonal )
{
    {
        ReadAccess<OtherType> rDiagonal( diagonal );
        WriteAccess<ValueType> wValues( mValues );
        ReadAccess<IndexType> rOffset( mOffset );

        for( IndexType i = 0; i < mNumRows; i++ )
        {
            for( IndexType ii = 0; ii < mNumDiagonals; ++ii )
            {
                const IndexType j = i + rOffset[ii];

                if( j >= mNumColumns )
                {
                    break;
                }

                if( j >= 0 )
                {
                    wValues[ii * mNumRows + i] *= static_cast<ValueType>( rDiagonal[j] );
                }
            }
        }
    }

    if( SCAI_LOG_TRACE_ON( logger ) )
    {
        SCAI_LOG_TRACE( logger, "DIA after scale diagonal" )
        print();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool DIAStorage<ValueType>::checkDiagonalProperty() const
{
    bool diagonalProperty = true;

    if( mNumRows != mNumColumns )
    {
        diagonalProperty = false;
    }
    else if( mNumRows == 0 )
    {
        // zero sized matrix has diagonal property

        diagonalProperty = true;
    }
    else if( mOffset.size() == 0 )
    {
        // full zero matrix but not zero size -> no diagonal property

        diagonalProperty = false;
    }
    else
    {
        // diagonal property is given if first diagonal is the main one

        ReadAccess<IndexType> offset( mOffset );
        diagonalProperty = offset[0] == 0;
    }

    SCAI_LOG_INFO( logger, *this << ": checkDiagonalProperty -> " << diagonalProperty )

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::check( const char* /* msg */) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumDiagonals, mOffset.size() )

    if( mNumDiagonals == 0 )
    {
        SCAI_ASSERT_EQUAL_ERROR( false, mDiagonalProperty )
    }

    SCAI_ASSERT_EQUAL_ERROR( mNumDiagonals * mNumRows, mValues.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_DEBUG( logger, "set identity, size = " << size )

    mNumRows = size;
    mNumColumns = size;

    mNumDiagonals = 1; // identity has exactly one diagonal

    ContextPtr loc = getContextPtr();

    {
        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )

        WriteOnlyAccess<IndexType> wOffset( mOffset, loc, mNumDiagonals );

        SCAI_CONTEXT_ACCESS( loc )

        setVal( wOffset.get(), 1, 0 );
    }

    {
        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )
        WriteOnlyAccess<ValueType> values( mValues, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        setVal( values.get(), mNumRows, static_cast<ValueType>(1.0) );
    }

    mDiagonalProperty = true;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setUsedDiagonal(
    bool upperDiagonalUsed[],
    bool lowerDiagonalUsed[],
    IndexType i,
    IndexType j )
{
    if( j >= i )
    {
        bool& flag = upperDiagonalUsed[j - i];

        // set flag only if not already set, improves cache usage

        if( !flag )
        {
            flag = true; // write only if not true,
        }
    }
    else
    {
        bool& flag = lowerDiagonalUsed[i - j];

        if( !flag )
        {
            flag = true; // write only if not true,
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::buildCSR(
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr /* loc */) const
{
    SCAI_REGION( "Storage.DIA->CSR" )

    // TODO all done on host, so loc is unused

    SCAI_LOG_INFO( logger,
                   "buildTypedCSRData<" << common::getScalarType<OtherValueType>() << ">" << " from DIA<" << common::getScalarType<ValueType>() << "> = " << *this << ", diagonal property = " << mDiagonalProperty )

    ReadAccess<IndexType> diaOffsets( mOffset );
    ReadAccess<ValueType> diaValues( mValues );

    WriteOnlyAccess<IndexType> csrIA( ia, mNumRows + 1 );

    // TODO: Check if eps = 0.0 is correct
    OpenMPDIAUtils::getCSRSizes( csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(),
                                 diaValues.get(), static_cast<ValueType>(0.0) );

    if( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = OpenMPCSRUtils::sizes2offsets( csrIA.get(), mNumRows );

    SCAI_LOG_INFO( logger, "CSR: #non-zero values = " << numValues )

    WriteOnlyAccess<IndexType> csrJA( *ja, numValues );
    WriteOnlyAccess<OtherValueType> csrValues( *values, numValues );

    // TOOD: Check if eps = 0.0 is correct
    OpenMPDIAUtils::getCSRValues( csrJA.get(), csrValues.get(), csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns,
                                  mNumDiagonals, diaOffsets.get(), diaValues.get(), static_cast<ValueType>(0.0) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType UNUSED( numValues ),
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    ContextPtr UNUSED( loc ) )
{
    SCAI_REGION( "Storage.DIA<-CSR" )

    // loc is ignored, we do it on the Host

    ReadAccess<IndexType> csrIA( ia );
    ReadAccess<IndexType> csrJA( ja );
    ReadAccess<OtherValueType> csrValues( values );

    _MatrixStorage::setDimension( numRows, numColumns );

    SCAI_LOG_DEBUG( logger, "fill DIA sparse matrix " << mNumRows << " x " << mNumColumns << " from csr data" )

    // build a set of all used lower and upper diagonals

    IndexType maxNumDiagonals = std::max( mNumRows, mNumColumns );

    scoped_array<bool> upperDiagonalUsed( new bool[maxNumDiagonals] );
    scoped_array<bool> lowerDiagonalUsed( new bool[maxNumDiagonals] );

    for( IndexType i = 0; i < maxNumDiagonals; i++ )
    {
        upperDiagonalUsed[i] = false;
        lowerDiagonalUsed[i] = false;
    }

    #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

    for( IndexType i = 0; i < mNumRows; ++i )
    {
        for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; jj++ )
        {
            IndexType j = csrJA[jj]; // column used
            setUsedDiagonal( upperDiagonalUsed.get(), lowerDiagonalUsed.get(), i, j );
        }
    }

    mDiagonalProperty = OpenMPCSRUtils::hasDiagonalProperty( numRows, csrIA.get(), csrJA.get() );

    // mDiagonalProperty forces upper diagonal to be the first one

    setOffsets( maxNumDiagonals, upperDiagonalUsed.get(), lowerDiagonalUsed.get() );

    // now we can allocate and set the values
    {
        ReadAccess<IndexType> offset( mOffset );
        WriteOnlyAccess<ValueType> myValues( mValues, mNumDiagonals * mNumRows );

        #pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)

        for( int i = 0; i < mNumRows; i++ )
        {
            for( IndexType d = 0; d < mNumDiagonals; d++ )
            {
                IndexType j = i + offset[d];

                if( j < 0 || j >= mNumColumns )
                {
                    continue;
                }

                ValueType& addrValue = myValues[diaindex( i, d, mNumRows, mNumDiagonals )];

                // check for j >= 0 and j < mNumColumns not needed here

                addrValue = static_cast<ValueType>(0.0);

                for( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    if( csrJA[jj] == j )
                    {
                        addrValue = static_cast<ValueType>( csrValues[jj] );
                        break;
                    }
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setOffsets(
    const IndexType maxNumDiagonals,
    const bool upperDiagonalUsed[],
    const bool lowerDiagonalUsed[] )
{
    // Help routine to set offsets from used diagonals

    mNumDiagonals = 0;

    IndexType firstIndex = 0;

    if( mDiagonalProperty )
    {
        firstIndex = 1;
        mNumDiagonals = 1;
    }

    for( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
    {
        if( upperDiagonalUsed[i] )
        {
            mNumDiagonals++;
        }

        if( lowerDiagonalUsed[i] )
        {
            mNumDiagonals++;
        }
    }

    SCAI_LOG_INFO( logger, "storage data requires " << mNumDiagonals << " diagonals a " << mNumRows << " values" )

    WriteOnlyAccess<IndexType> wOffset( mOffset, mNumDiagonals );

    if( mNumDiagonals > 0 )
    {
        mNumDiagonals = 0;
        firstIndex = 0;

        if( mDiagonalProperty )
        {
            wOffset[mNumDiagonals++] = 0;
            firstIndex = 1;
        }

        for( IndexType i = maxNumDiagonals - 1; i > 0; i-- )
        {
            if( lowerDiagonalUsed[i] )
            {
                wOffset[mNumDiagonals++] = -i;
            }
        }

        SCAI_LOG_INFO( logger, "lower diagonals = " << mNumDiagonals )

        for( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
        {
            if( upperDiagonalUsed[i] )
            {
                wOffset[mNumDiagonals++] = i;
            }
        }

        SCAI_LOG_INFO( logger, "lower + upper diagonals = " << mNumDiagonals )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumDiagonals, wOffset.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::~DIAStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~DIAStorage for matrix " << mNumRows << " x " << mNumColumns << ", # diags = " << mNumDiagonals )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::purge()
{
    mNumColumns = 0;
    mNumRows = 0;
    mNumDiagonals = 0;

    mOffset.purge();
    mValues.purge();

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate DIA sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    mNumRows = numRows;
    mNumColumns = numColumns;

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DIAStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << mNumRows << " x " << mNumColumns
           << ", nd = " << mNumDiagonals << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( asum, loc, BLAS, BLAS1, ValueType )

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return asum( mValues.size(), data.get(), 1, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( dot, loc, BLAS, BLAS1, ValueType )

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return ::sqrt(dot( mValues.size(), data.get(), 1, data.get(), 1, NULL ));
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, DIAUtils, Reductions, ValueType )

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal( mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(), diaValues.get() );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_LOG_DEBUG( logger, "get value (" << i << ", " << j << ") from " << *this )

    const ReadAccess<IndexType> offset( mOffset );

    ValueType myValue = static_cast<ValueType>(0.0);

    // check for a matching diagonal element in the row i

    for( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        if( i + offset[d] == j )
        {
            const ReadAccess<ValueType> values( mValues );
            SCAI_LOG_DEBUG( logger,
                            "get value (" << i << ", " << j << ") is diag = " << d << ", offset = " << offset[d] << ", index = " << diaindex( i, d, mNumRows, mNumColumns ) )
            myValue = values[diaindex( i, d, mNumRows, mNumDiagonals )];
            break;
        }
    }

    return myValue;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mOffset.prefetch( location );
    mValues.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DIAStorage<ValueType>::getNumDiagonals() const
{
    return mNumDiagonals;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& DIAStorage<ValueType>::getOffsets() const
{
    return mOffset;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<ValueType>& DIAStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::wait() const
{
    mOffset.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::swap( DIAStorage<ValueType>& other )
{
    std::swap( mNumDiagonals, other.mNumDiagonals );
    mOffset.swap( other.mOffset );
    mValues.swap( other.mValues );

    MatrixStorage<ValueType>::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t DIAStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType) * mOffset.size();
    memoryUsage += sizeof(ValueType) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::matrixTimesVector(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.DIA.timesVector" )

    ContextPtr loc = getContextPtr();

    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result << " on " << *loc )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    LAMA_INTERFACE_FN_DEFAULT_T( normalGEMV, loc, DIAUtils, Mult, ValueType )

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y is handled by coressponding accesses

    if( &result == &y )
    {
        SCAI_LOG_DEBUG( logger, "result == y" )

        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets.get(), diaValues.get(), NULL );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "result != y" )

        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets.get(), diaValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::vectorTimesMatrix(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.DIA.VectorTimesMatrix" )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    ContextPtr loc = getContextPtr();

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )

    LAMA_INTERFACE_FN_DEFAULT_T( normalGEVM, loc, DIAUtils, Mult, ValueType )

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEVM can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets.get(), diaValues.get(), NULL );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets.get(), diaValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::matrixTimesVectorAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.DIA.timesVectorAsync" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( normalGEMV, loc, DIAUtils, Mult, ValueType )

    if( loc->getType() == context::Host )
    {
        // Start directly a task, avoids pushing of accesses

        void (DIAStorage::*mv)(
            LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>& ) const

            = &DIAStorage<ValueType>::matrixTimesVector;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        return new TaskSyncToken( bind( mv, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    // logging + checks not needed when started as a task

    SCAI_LOG_INFO( logger,
                   "Start z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result << " on " << *loc )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > diaOffsets( new ReadAccess<IndexType>( mOffset, loc ) );
    shared_ptr<ReadAccess<ValueType> > diaValues( new ReadAccess<ValueType>( mValues, loc ) );

    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y is handled by coressponding accesses

    if( &result == &y )
    {
        SCAI_LOG_DEBUG( logger, "result == y" )

        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        syncToken->pushToken( wResult );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets->get(), diaValues->get(), syncToken.get() );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "result != y" )

        shared_ptr<WriteOnlyAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumRows ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        syncToken->pushToken( rY );
        syncToken->pushToken( wResult );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets->get(), diaValues->get(), syncToken.get() );
    }

    syncToken->pushToken( rX );
    syncToken->pushToken( diaValues );
    syncToken->pushToken( diaOffsets );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::vectorTimesMatrixAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.DIA.vectorTimesMatrixAsync" )

    ContextPtr loc = getContextPtr();

    // Note: checks will be done by asynchronous task in any case
    //       and exception in tasks are handled correctly

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on " << *loc )

    if( loc->getType() == context::Host )
    {
        // execution as separate thread

        void (DIAStorage::*pf)(
            LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>& ) const

            = &DIAStorage<ValueType>::vectorTimesMatrix;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on Host by own thread" )

        return new TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    LAMA_INTERFACE_FN_T( normalGEVM, loc, DIAUtils, Mult, ValueType )

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > diaOffsets( new ReadAccess<IndexType>( mOffset, loc ) );
    shared_ptr<ReadAccess<ValueType> > diaValues( new ReadAccess<ValueType>( mValues, loc ) );

    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets->get(), diaValues->get(), syncToken.get() );

        syncToken->pushToken( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumColumns ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumColumns, mNumDiagonals,
                    diaOffsets->get(), diaValues->get(), syncToken.get() );

        syncToken->pushToken( wResult );
        syncToken->pushToken( rY );
    }

    syncToken->pushToken( rX );
    syncToken->pushToken( diaValues );
    syncToken->pushToken( diaOffsets );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
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

    WriteAccess<ValueType> wSolution( solution );
    ReadAccess<IndexType> diaOffset( mOffset );
    ReadAccess<ValueType> diaValues( mValues );
    ReadAccess<ValueType> rOldSolution( oldSolution );
    ReadAccess<ValueType> rRhs( rhs );

    OpenMPDIAUtils::jacobi( wSolution.get(), mNumColumns, mNumDiagonals, diaOffset.get(), diaValues.get(),
                            rOldSolution.get(), rRhs.get(), omega, mNumRows, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::clone() const
{
    return new DIAStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::copy() const
{
    return new DIAStorage<ValueType>( *this );
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

#define LAMA_DIA_STORAGE_INSTANTIATE(z, I, _)                                     \
    template<>                                                                    \
    const char* DIAStorage<ARITHMETIC_HOST_TYPE_##I>::typeName()                  \
    {                                                                             \
        return "DIAStorage<" PRINT_STRING(ARITHMETIC_HOST_TYPE_##I) ">";      \
    }                                                                             \
                                                                                  \
    template class COMMON_DLL_IMPORTEXPORT DIAStorage<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_DIA_STORAGE_INSTANTIATE, _ )

#undef LAMA_DIA_STORAGE_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
