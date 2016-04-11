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

// internal scai libraries
#include <scai/sparsekernel/DIAKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::scoped_array;
using common::shared_ptr;

using tasking::SyncToken;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;

using sparsekernel::CSRKernelTrait;
using sparsekernel::DIAKernelTrait;

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
    const HArray<IndexType>& offsets,
    const HArray<ValueType>& values )

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
Format::MatrixStorageFormat DIAStorage<ValueType>::getFormat() const
{
    return Format::DIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": has not diagonal property, cannot set diagonal" )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;

    // take context of this storage to set

    ContextPtr loc = setVal.getValidContext( this->getContextPtr() ); 

    {
        // not all values might be changed, so use WriteAccess instead of WriteOnlyAccess

        WriteAccess<ValueType> wValues( mValues, loc );
        ReadAccess<IndexType> rOffset( mOffset, loc );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[loc]( wValues.get(), numDiagonalElements, value, common::reduction::COPY );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType i ) const
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
void DIAStorage<ValueType>::getDiagonalImpl( HArray<OtherType>& diagonal ) const
{
    static LAMAKernel<UtilKernelTrait::set<OtherType, ValueType> > set;

    ContextPtr loc = set.getValidContext( getContextPtr() );

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // Diagonal is first column

    set[ loc ]( wDiagonal.get(), rValues.get(), numDiagonalElements, common::reduction::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::setDiagonalImpl( const HArray<OtherType>& diagonal )
{
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    numDiagonalElements = std::min( numDiagonalElements, diagonal.size() );

    static LAMAKernel<UtilKernelTrait::set<ValueType, OtherType> > set;

    ContextPtr loc = set.getValidContext( getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    set[loc]( wValues.get(), rDiagonal.get(), numDiagonalElements, common::reduction::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scaleImpl( const ValueType value )
{
    HArrayUtils::scale( mValues, value, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::conj()
{
    HArrayUtils::conj( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::scaleImpl( const HArray<OtherType>& diagonal )
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

    {
        static LAMAKernel<UtilKernelTrait::setVal<IndexType, IndexType> > setVal;

        ContextPtr loc = setVal.getValidContext( getContextPtr() );

        WriteOnlyAccess<IndexType> wOffset( mOffset, loc, mNumDiagonals );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[ loc ]( wOffset.get(), 1, 0, common::reduction::COPY );
    }

    {
        static LAMAKernel<UtilKernelTrait::setVal<ValueType, ValueType> > setVal;

        ContextPtr loc = setVal.getValidContext( getContextPtr() );

        WriteOnlyAccess<ValueType> values( mValues, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        setVal[ loc ]( values.get(), mNumRows, ValueType( 1 ), common::reduction::COPY );
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
template<typename CSRValueType>
void DIAStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<CSRValueType>* values,
    const ContextPtr prefLoc ) const
{
    SCAI_REGION( "Storage.DIA->CSR" )

    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;
    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, CSRValueType> > getCSRValues;

    // do it where all routines are avaialble

    ContextPtr loc = sizes2offsets.getValidContext( getCSRSizes, getCSRValues, prefLoc );

    SCAI_LOG_INFO( logger,
                   "buildTypedCSRData<" << common::getScalarType<CSRValueType>() << ">" 
                    << " from DIA<" << common::getScalarType<ValueType>() << "> = " << *this << ", diagonal property = " << mDiagonalProperty )

    ReadAccess<IndexType> diaOffsets( mOffset );
    ReadAccess<ValueType> diaValues( mValues );

    WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );

    // In contrary to COO and CSR, the DIA format stores also some ZERO values like Dense

    ValueType eps = static_cast<ValueType>( 0.0 );

    getCSRSizes[loc]( csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(),
                      diaValues.get(), eps );

    if( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets[loc]( csrIA.get(), mNumRows );

    SCAI_LOG_INFO( logger, "CSR: #non-zero values = " << numValues )

    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<CSRValueType> csrValues( *values, loc, numValues );

    getCSRValues[loc]( csrJA.get(), csrValues.get(), csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns,
                       mNumDiagonals, diaOffsets.get(), diaValues.get(), eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType UNUSED( numValues ),
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    ContextPtr UNUSED( prefLoc ) )
{
    SCAI_REGION( "Storage.DIA<-CSR" )

    static LAMAKernel<CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;

    // prefLoc is ignored, we do it on the Host
    // ToDo: replace Host code with kernels, implement kernels for other devices

    ContextPtr loc = Context::getHostPtr();

    ReadAccess<IndexType> csrIA( ia, loc );
    ReadAccess<IndexType> csrJA( ja, loc );
    ReadAccess<OtherValueType> csrValues( values, loc );

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

    mDiagonalProperty = hasDiagonalProperty[loc]( numRows, csrIA.get(), csrJA.get() );

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

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;

    ContextPtr loc = asum.getValidContext( this->getContextPtr() );

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return asum[loc]( mValues.size(), data.get(), 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;

    ContextPtr loc = dot.getValidContext( this->getContextPtr() );

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return common::Math::sqrt(dot[loc]( mValues.size(), data.get(), 1, data.get(), 1 ));
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    static LAMAKernel<DIAKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = absMaxVal.getValidContext( this->getContextPtr() );

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal[loc]( mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(), diaValues.get() );

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
const HArray<IndexType>& DIAStorage<ValueType>::getOffsets() const
{
    return mOffset;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& DIAStorage<ValueType>::getValues() const
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
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.DIA.timesVector" )

    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y" 
                    << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )

    if ( beta != common::constants::ZERO )
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
    }

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = normalGEMV.getValidContext( this->getContextPtr() );

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    // Note: read access to y must appear before write access to result in case of alias

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );  // result might be aliased to y

    SCAI_CONTEXT_ACCESS( loc )

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                     diaOffsets.get(), diaValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.DIA.VectorTimesMatrix" )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    static LAMAKernel<DIAKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = normalGEVM.getValidContext( this->getContextPtr() );

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

    SCAI_CONTEXT_ACCESS( loc )

    normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                     diaOffsets.get(), diaValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.DIA.timesVectorAsync" )

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = normalGEMV.getValidContext( this->getContextPtr() );

    // logging + checks not needed when started as a task

    SCAI_LOG_INFO( logger,
                   "Start z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result << " on " << *loc )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    SCAI_ASYNCHRONOUS( *syncToken ) 

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

    SCAI_CONTEXT_ACCESS( loc )

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                     diaOffsets.get(), diaValues.get() );

    syncToken->pushRoutine( rY.releaseDelayed() );
    syncToken->pushRoutine( wResult.releaseDelayed() );
    syncToken->pushRoutine( rX.releaseDelayed() );
    syncToken->pushRoutine( diaValues.releaseDelayed() );
    syncToken->pushRoutine( diaOffsets.releaseDelayed() );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::vectorTimesMatrixAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.DIA.vectorTimesMatrixAsync" )

    static LAMAKernel<DIAKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = normalGEVM.getValidContext( this->getContextPtr() );

    // Note: checks will be done by asynchronous task in any case
    //       and exception in tasks are handled correctly

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on " << *loc )

    if( loc->getType() == Context::Host )
    {
        // execution as separate thread

        void (DIAStorage::*pf)(
            HArray<ValueType>&,
            const ValueType,
            const HArray<ValueType>&,
            const ValueType,
            const HArray<ValueType>& ) const

            = &DIAStorage<ValueType>::vectorTimesMatrix;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on Host by own thread" )

        return new tasking::TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    syncToken->setCurrent();

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues(  mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if( &result == &y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumDiagonals,
                         diaOffsets.get(), diaValues.get() );

        syncToken->pushRoutine( wResult.releaseDelayed() );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                         diaOffsets.get(), diaValues.get() );

        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
    }

    syncToken->pushRoutine( rX.releaseDelayed() );
    syncToken->pushRoutine( diaValues.releaseDelayed() );
    syncToken->pushRoutine( diaOffsets.releaseDelayed() );

    syncToken->unsetCurrent();

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )

    // matrix must be square

    static LAMAKernel<DIAKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = jacobi.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> diaOffset( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );

    jacobi[loc]( wSolution.get(), mNumColumns, mNumDiagonals, diaOffset.get(), diaValues.get(),
                 rOldSolution.get(), rRhs.get(), omega, mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::copy() const
{
    return new DIAStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::newMatrixStorage() const
{
   common::unique_ptr<DIAStorage<ValueType> > storage( new DIAStorage<ValueType>() ); 
   storage->setContextPtr( this->getContextPtr() );
   return storage.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* DIAStorage<ValueType>::create()
{
    return new DIAStorage<ValueType>();
}

template<typename ValueType>
std::string DIAStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string("DIAStorage<") << common::getScalarType<ValueType>() << std::string(">");
    return s.str();
}

template<typename ValueType>
const char* DIAStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType DIAStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DIA, common::getScalarType<ValueType>() );
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DIAStorage, ARITHMETIC_HOST_CNT, ARITHMETIC_HOST )

} /* end namespace lama */

} /* end namespace scai */
