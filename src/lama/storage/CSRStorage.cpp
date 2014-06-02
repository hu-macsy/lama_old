/**
 * @file CSRStorage.cpp
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
 * @brief Implementation and instantiation for template class CSRStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 * @since 1.0.0
 */

// hpp
#include <lama/storage/CSRStorage.hpp>

// others
#include <lama/LAMAInterface.hpp>
#include <lama/ContextAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/HostWriteAccess.hpp>

#include <lama/storage/StorageMethods.hpp>

#include <lama/LAMAArrayUtils.hpp>

#include <lama/distribution/Redistributor.hpp>

#include <lama/task/TaskSyncToken.hpp>

#include <lama/openmp/OpenMPUtils.hpp>
#include <lama/openmp/OpenMPBLAS1.hpp>
#include <lama/openmp/OpenMPCSRUtils.hpp>

// assert
#include <lama/exception/LAMAAssert.hpp>

// boost
#include <boost/bind.hpp>
#include <boost/preprocessor.hpp>
#include <lama/tracing.hpp>

#include <cmath>

namespace lama
{

using std::abs;
using std::auto_ptr;
using boost::shared_ptr;

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, CSRStorage<ValueType>::logger, "MatrixStorage.CSRStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage() :
    CRTPMatrixStorage<CSRStorage<ValueType>, ValueType> (0, 0),
    mNumValues( 0 ),
    mSortedRows( false )
{
    allocate( 0, 0 ); // creates at least mIa
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const _LAMAArray& values )

    : CRTPMatrixStorage<CSRStorage<ValueType>,ValueType>()
{
    this->setCSRData( numRows, numColumns, numValues, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::print() const
{
    using std::cout;
    using std::endl;

    cout << "CSRStorage " << mNumRows << " x " << mNumColumns << ", #values = " << mNumValues << endl;

    HostReadAccess<IndexType> ia( mIa );
    HostReadAccess<IndexType> ja( mJa );
    HostReadAccess<ValueType> values( mValues );

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        cout << "Row " << i << " ( " << ia[i] << " - " << ia[i + 1] << " ) :";

        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            cout << " " << ja[jj] << ":" << values[jj];
        }
        cout << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;

    mIa.clear();
    mJa.clear();
    mValues.clear();

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat CSRStorage<ValueType>::getFormat() const
{
    return CSR;
}

/* --------------------------------------------------------------------------- */

#ifndef LAMA_ASSERT_LEVEL_DEBUG
template<typename ValueType>
void CSRStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void CSRStorage<ValueType>::check( const char* msg ) const
{
    LAMA_ASSERT_EQUAL_ERROR( mNumRows + 1, mIa.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mJa.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )

    // check ascending values in offset array mIa

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT_T( isSorted, loc, Utils, Reductions, IndexType )
        LAMA_INTERFACE_FN_T( getValue, loc, Utils, Getter, IndexType )

        ReadAccess<IndexType> csrIA( mIa, loc );

        LAMA_CONTEXT_ACCESS( loc )

        bool ascending = true;  // check for ascending

        IndexType numValues = getValue( csrIA.get(), mNumRows );
 
        LAMA_ASSERT_ERROR( numValues == mNumValues, "ia[" << mNumRows << "] = " << numValues <<
                             ", expected " << mNumValues << ", msg = " << msg )

        LAMA_ASSERT_ERROR( isSorted ( csrIA.get(), mNumRows + 1, ascending ),
                           *this << " @ " << msg << ": IA is illegal offset array" )
    }

    // check column indexes in JA

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT( validIndexes, loc, Utils, Indexes )

        ReadAccess<IndexType> rJA( mJa, loc );

        LAMA_CONTEXT_ACCESS( loc )

        LAMA_ASSERT_ERROR( validIndexes ( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool CSRStorage<ValueType>::checkDiagonalProperty() const
{
    // diagonal property is given if size of matrix is 0

    if ( mNumRows == 0  || mNumColumns == 0 )
    {
        return true;
    }

    // non-zero sized matrix with no values has not diagonal property

    if ( mNumValues == 0 )
    {
        return false;
    }

    ContextPtr loc = getContextPtr();  // there we do the checks

    //get function pointer
    LAMA_INTERFACE_FN( hasDiagonalProperty, loc, CSRUtils, Offsets )

    //get read access
    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    LAMA_CONTEXT_ACCESS( loc )

    IndexType numDiagonals = std::min( mNumRows, mNumColumns );

    bool diagonalProperty = hasDiagonalProperty( numDiagonals, csrIA.get(), csrJA.get() );

    LAMA_LOG_DEBUG( logger, *this << ": diagonalProperty = " << diagonalProperty );

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setIdentity( const IndexType size )
{
    LAMA_LOG_DEBUG( logger, "set identity, size = " << size )

    _MatrixStorage::setDimension( size, size );

    mNumValues = mNumRows;

    HostWriteOnlyAccess<IndexType> ia( mIa, mNumRows + 1 );
    HostWriteOnlyAccess<IndexType> ja( mJa, mNumValues );
    HostWriteOnlyAccess<ValueType> values( mValues, mNumValues );

    ValueType one = static_cast<ValueType>( 1.0 );

    OpenMPUtils::setOrder( ia.get(), mNumRows + 1 );
    OpenMPUtils::setOrder( ja.get(), mNumRows );
    OpenMPUtils::setVal( values.get(), mNumRows, one );

    mDiagonalProperty = true; // obviously given for identity matrix
    mSortedRows = true; // obviously given for identity matrix

    // Note: we do not build row indexes, no row is empty

    LAMA_LOG_INFO( logger, *this << ": identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    const ContextPtr /* loc */)
{
    ContextPtr loc = this->getContextPtr();

    if ( ia.size() == numRows )
    {
        // checking is done where ia is already valid, preferred is loc

        ContextPtr loc1 = ia.getValidContext( loc->getType() );

        // we assume that ia contains the sizes, verify it by summing up

        LAMA_INTERFACE_FN_DEFAULT_T( sum, loc1, Utils, Reductions, IndexType )

        ReadAccess<IndexType> csrIA( ia, loc1 );
        
        LAMA_CONTEXT_ACCESS( loc1 )
        
        IndexType n = sum( csrIA.get(), numRows );
 
        if ( n != numValues )
        {
            LAMA_THROWEXCEPTION( "ia is invalid size array" )
        }
    }
    else if ( ia.size() == numRows + 1 )
    {
        // checking is done where ia is already valid

        ContextPtr loc1 = ia.getValidContext( loc->getType() );

        LAMA_INTERFACE_FN_DEFAULT( validOffsets, loc1, CSRUtils, Offsets )

        ReadAccess<IndexType> csrIA( ia, loc1 );
    
        LAMA_CONTEXT_ACCESS( loc1 )
        
        if ( !validOffsets( csrIA.get(), numRows, numValues ) )
        {
            LAMA_THROWEXCEPTION( "ia is invalid offset array" )
        }
    }
    else
    {
        LAMA_THROWEXCEPTION( "ia array with size = " << ia.size() << " illegal, #rows = " << numRows )
    }

    LAMA_ASSERT_EQUAL_ERROR( numValues, ja.size() );
    LAMA_ASSERT_EQUAL_ERROR( numValues, values.size() );

    {
        LAMA_INTERFACE_FN( validIndexes, loc, Utils, Indexes )

        // make sure that column indexes in JA are all valid 

        ReadAccess<IndexType> csrJA( ja, loc );

        LAMA_CONTEXT_ACCESS( loc )

        if ( !validIndexes( csrJA.get(), numValues, numColumns ) )
        {
            LAMA_THROWEXCEPTION( "invalid column indexes in ja = " << ja << ", #columns = " << numColumns )
        }
    }

    // now we can copy all data

    mNumRows = numRows;
    mNumColumns = numColumns;
    mNumValues = numValues;

    LAMA_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    // storage data will be directly allocated on the location

    if ( ia.size() == numRows )
    {
        {
            // reserve enough memory for mIa

            WriteOnlyAccess<IndexType> myIA( mIa, loc, mNumRows + 1 );
        }

        LAMAArrayUtils::assign( mIa, ia, loc );
            
        {
            ContextPtr loc1 = loc; // loc might change if sizes2offsets is not available

            LAMA_INTERFACE_FN_DEFAULT( sizes2offsets, loc1, CSRUtils, Offsets )

            WriteAccess<IndexType> myIA( mIa, loc1 );

            myIA.resize( mNumRows + 1 );                  // no realloc as capacity is sufficient

            sizes2offsets( myIA.get(), numRows );
        }
    }
    else
    {
        LAMAArrayUtils::assign( mIa, ia, loc );
    }

    LAMAArrayUtils::assign( mValues, values, loc );
    LAMAArrayUtils::assign( mJa, ja, loc );

    /* do not sort rows, destroys diagonal property during redistribute

     OpenMPCSRUtils::sortRowElements( myJA.get(), myValues.get(), myIA.get(),
     mNumRows, mDiagonalProperty );

     */

    mDiagonalProperty = checkDiagonalProperty();

    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::sortRows( bool diagonalProperty )
{
    {
        HostReadAccess<IndexType> csrIA( mIa );
        HostWriteAccess<IndexType> csrJA( mJa );
        HostWriteAccess<ValueType> csrValues( mValues );

        OpenMPCSRUtils::sortRowElements( csrJA.get(), csrValues.get(), csrIA.get(),
                                         mNumRows, diagonalProperty );
    }

    mDiagonalProperty = checkDiagonalProperty();
}

//this version avoids copying the ia, ja, and value arrays, but instead swaps them
//also it does not check their validity
//this is much faster of course, but destroys the input ia, ja and value arrays
template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::setCSRDataSwap(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>& ja,
    LAMAArray<OtherValueType>& values,
    const ContextPtr /* loc */ )
{
    //set necessary information
    mNumRows = numRows;
    mNumColumns = numColumns;
    mNumValues = numValues;

    LAMA_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    //swap arrays

    mIa.swap( ia );
    mJa.swap( ja );
    mValues.swap( values );

    mDiagonalProperty = checkDiagonalProperty();

    // this builds only row indices if context is on host

    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildRowIndexes()
{
    LAMA_LOG_INFO( logger, "buildRowIndexes is temporarily disabled." );
    return;

    mRowIndexes.clear();

    if ( mNumRows == 0 )
    {
        return;
    }

    if ( getContext().getType() != Context::Host )
    {
        LAMA_LOG_INFO( logger, "CSRStorage: build row indices is currently only implemented on host" )
    }

    // This routine is only available on the Host

    ContextPtr loc = ContextFactory::getContext( Context::Host );

    ReadAccess<IndexType> csrIA( mIa, loc );

    IndexType nonZeroRows = OpenMPCSRUtils::countNonEmptyRowsByOffsets( csrIA.get(), mNumRows );

    float usage = float( nonZeroRows ) / float( mNumRows );

    if ( usage >= mCompressThreshold )
    {
        LAMA_LOG_INFO( logger, "CSRStorage: do not build row indexes, usage = " << usage )
        return;
    }

    LAMA_LOG_INFO( logger, "CSRStorage: build row indexes, #entries = " << nonZeroRows )

    WriteOnlyAccess<IndexType> rowIndexes( mRowIndexes, loc, nonZeroRows );

    OpenMPCSRUtils::setNonEmptyRowsByOffsets( rowIndexes.get(), nonZeroRows, csrIA.get(), mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::redistributeCSR( const CSRStorage<ValueType>& other, const Redistributor& redistributor )
{
    LAMA_REGION( "Storage.redistributeCSR" )

    const Distribution& sourceDistribution = *redistributor.getSourceDistributionPtr();
    const Distribution& targetDistribution = *redistributor.getTargetDistributionPtr();

    LAMA_LOG_INFO( logger, other << ": redistribute of CSR<" << other.getValueType()
                           << "> to CSR<" << this->getValueType() << " via " << redistributor )

    bool sameDist = false;

    // check for same distribution, either equal or both replicated

    if ( sourceDistribution.isReplicated() && targetDistribution.isReplicated() )
    {
        sameDist = true;
    }
    else if ( &sourceDistribution == &targetDistribution )
    {
        sameDist = true;
    }

    if ( sameDist )
    {
        LAMA_LOG_INFO( logger, "redistributor with same source/target distribution" )

        assign( other );

        return; // so we are done
    }

    // check that source distribution fits with storage

    LAMA_ASSERT_EQUAL_ERROR( other.getNumRows(), sourceDistribution.getLocalSize() )

    if ( &other == this )
    {
        // due to alias we need temporary array

        LAMAArray<IndexType> targetIA;
        LAMAArray<IndexType> targetJA;
        LAMAArray<ValueType> targetValues;

        StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues,
                                                    other.getIA(), other.getJA(), other.getValues(),
                                                    redistributor );

        // we can swap the new arrays

        mIa.swap( targetIA );
        mJa.swap( targetJA );
        mValues.swap( targetValues );
    }
    else
    {
        StorageMethods<ValueType>::redistributeCSR( mIa, mJa, mValues,
                                                    other.getIA(), other.getJA(), other.getValues(),
                                                    redistributor );
    }

    // it is not necessary to convert the other storage to CSR


    mNumColumns = other.getNumColumns();
    mNumRows = mIa.size() - 1;
    mNumValues = mJa.size();

    mDiagonalProperty = checkDiagonalProperty();

    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::~CSRStorage()
{
    LAMA_LOG_DEBUG( logger,
                    "~CSRStorage for maxtrix " << mNumRows << " x " << mNumColumns << ", # non-zeros = " << mNumValues )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType CSRStorage<ValueType>::getNumValues() const
{
    return mNumValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::purge()
{
    mNumColumns = 0;
    mNumRows = 0;
    mNumValues = 0;

    mIa.purge();
    mJa.purge();
    mValues.purge();

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    LAMA_LOG_INFO( logger,
                   "allocate CSR sparse matrix of size " << numRows << " x " << numColumns << ", numValues = 0" )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValues = 0;

    mJa.clear();
    mValues.clear();

    HostWriteOnlyAccess<IndexType> ia( mIa, mNumRows + 1 );

    // make a correct initialization for the offset array

    OpenMPUtils::setVal( ia.get(), mNumRows + 1, 0 );

    mDiagonalProperty = false;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::compress( const ValueType eps /* = 0.0 */)
{
    HostWriteAccess<IndexType> ia( mIa );
    HostReadAccess<IndexType> ja( mJa );
    HostReadAccess<ValueType> values( mValues );

    IndexType nonDiagZeros = 0;

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            if ( ja[jj] == i )
            {
                continue;
            }

            if ( abs( values[jj] ) <= eps )
            {
                ++nonDiagZeros;
            }
        }
    }

    LAMA_LOG_INFO( logger, "compress: " << nonDiagZeros << " non-diagonal zero elements" )

    if ( nonDiagZeros == 0 )
    {
        return;
    }

    const IndexType newNumValues = mJa.size() - nonDiagZeros;

    LAMAArray<ValueType> newValuesArray;
    LAMAArray<IndexType> newJaArray;

    HostWriteOnlyAccess<ValueType> newValues( newValuesArray, newNumValues );
    HostWriteOnlyAccess<IndexType> newJa( newJaArray, newNumValues );

    IndexType gap = 0;

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        for ( IndexType jj = ia[i] + gap; jj < ia[i + 1]; ++jj )
        {
            if ( abs( values[jj] ) <= eps && ja[jj] != i )
            {
                ++gap;
                continue;
            }

            newValues[jj - gap] = values[jj];
            newJa[jj - gap] = ja[jj];
        }
        ia[i + 1] -= gap;
    }

    LAMA_ASSERT_EQUAL_DEBUG( gap, nonDiagZeros )

    ia.release();
    ja.release();
    values.release();
    newJa.release();
    newValues.release();

    mJa.swap( newJaArray );
    mValues.swap( newValuesArray );
    mNumValues = newNumValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::swap( CSRStorage<ValueType>& other )
{
    std::swap( mNumValues, other.mNumValues );
    mIa.swap( other.mIa );
    mJa.swap( other.mJa );
    mValues.swap( other.mValues );

    // swap sizes and row indexes

    MatrixStorage<ValueType>::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::swap( LAMAArray<IndexType>& ia, LAMAArray<IndexType>& ja, LAMAArray<ValueType>& values )
{
    LAMA_ASSERT_EQUAL_ERROR( ia.size(), mNumRows + 1 )

    IndexType numValues = 0;

    {
        HostReadAccess<IndexType> csrIA( ia );
        numValues = csrIA[mNumRows];
    }

    LAMA_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    LAMA_ASSERT_EQUAL_ERROR( numValues, values.size() )

    mNumValues = numValues;

    mIa.swap( ia );
    mJa.swap( ja );
    mValues.swap( values );

    mDiagonalProperty = checkDiagonalProperty();

    // build new array of row indexes

    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t CSRStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType) * mIa.size();
    memoryUsage += sizeof(IndexType) * mJa.size();
    memoryUsage += sizeof(ValueType) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "CSR( rows = " << mNumRows << ", cols= " << mNumColumns << ", nz = " << mNumValues << ", diag = "
           << mDiagonalProperty << ", sorted = " << mSortedRows << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    LAMA_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )
    const HostReadAccess<IndexType> ia( mIa );
    const HostReadAccess<IndexType> ja( mJa );
    const HostReadAccess<ValueType> values( mValues );
    ValueType myValue = 0;

    LAMA_LOG_TRACE( logger, "search column in ja from " << ia[i] << ":" << ia[i + 1] )

    for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
    {
        IndexType col = ja[jj];

        LAMA_ASSERT_DEBUG( 0 <= col && col < mNumColumns,
                           "column index at pos " << jj << " = " << col << " out of range" )

        if ( col == j )
        {
            LAMA_LOG_TRACE( logger, "found column j = " << j << " at " << jj << ", value = " << values[jj] )
            myValue = values[jj];
            break;
        }
    }

    return myValue;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mIa.prefetch( location );
    mJa.prefetch( location );
    mValues.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
LAMAArray<IndexType>& CSRStorage<ValueType>::getIA()
{
    return mIa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
LAMAArray<IndexType>& CSRStorage<ValueType>::getJA()
{
    return mJa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
LAMAArray<ValueType>& CSRStorage<ValueType>::getValues()
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& CSRStorage<ValueType>::getIA() const
{
    return mIa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& CSRStorage<ValueType>::getJA() const
{
    return mJa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<ValueType>& CSRStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setDiagonalImpl( const Scalar value )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    if ( !mDiagonalProperty )
    {
        LAMA_THROWEXCEPTION( "setDiagonal: matrix storage has not diagonal property." )
    }

    ValueType val = value.getValue<ValueType>();

    HostReadAccess<IndexType> wIa( mIa );
    HostWriteAccess<ValueType> wValues( mValues );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wValues[wIa[i]] = val;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherValueType>& diagonal )
{
    IndexType numDiagonalElements = diagonal.size();

    {
        HostReadAccess<OtherValueType> rDiagonal( diagonal );
        HostReadAccess<IndexType> csrIA( mIa );

        HostWriteAccess<ValueType> wValues( mValues ); // partial setting

        //  wValues[ wIa[ i ] ] = rDiagonal[ i ];

        OpenMPUtils::setScatter( wValues.get(), csrIA.get(), rDiagonal.get(), numDiagonalElements );
    }

    if ( LAMA_LOG_TRACE_ON( logger ) )
    {
        LAMA_LOG_TRACE( logger, "CSR after setDiagonal" )
        print();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void CSRStorage<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    LAMA_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    HostWriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const HostReadAccess<IndexType> ia( mIa );
    const HostReadAccess<IndexType> ja( mJa );
    const HostReadAccess<ValueType> values( mValues );

    for ( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = 0.0;
    }
    for ( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
    {
        wRow[ja[jj]] = static_cast<OtherType>( values[jj] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherValueType>& diagonal ) const
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
   
    ContextPtr loc = getContextPtr(); 
    WriteAccess<OtherValueType> wDiagonal( diagonal, loc, numDiagonalElements, false );
    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<ValueType> rValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )
    LAMA_INTERFACE_FN_TT( setGather, loc, Utils, Copy, OtherValueType, ValueType )

    setGather( wDiagonal.get(), rValues.get(), csrIA.get(), numDiagonalElements );

/*
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    HostWriteOnlyAccess<OtherValueType> wDiagonal( diagonal, numDiagonalElements );

    HostReadAccess<IndexType> csrIA( mIa );
    HostReadAccess<ValueType> rValues( mValues );

    //  diagonal[ i ] = rValues[ csrIA[ i ] ], diagonal values are at the offsets

    OpenMPUtils::setGather( wDiagonal.get(), rValues.get(), csrIA.get(), numDiagonalElements );
*/
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( scale, loc, Utils, Transform, ValueType )

    WriteAccess<ValueType> csrValues( mValues, loc );

    ValueType value = scalar.getValue<ValueType>();

    scale( csrValues.get(), value, mNumValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::scaleImpl( const LAMAArray<OtherValueType>& diagonal )
{
    IndexType n = std::min( mNumRows, diagonal.size() );
   
    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( scaleRows, loc, CSRUtils, Scale, ValueType, OtherValueType )

    LAMA_CONTEXT_ACCESS( loc )

    {
        ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
        ReadAccess<IndexType> csrIA( mIa, loc );
        WriteAccess<ValueType> csrValues( mValues, loc ); // updateAccess

        scaleRows( csrValues.get(), csrIA.get(), n, rDiagonal.get() );
    }

    if ( LAMA_LOG_TRACE_ON( logger ) )
    {
        LAMA_LOG_TRACE( logger, "CSR after scale diagonal" )
        print();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::wait() const
{
    mIa.wait();
    mJa.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::assign( const _MatrixStorage& other )
{
    if ( &other == this )
    {
        // this special case avoids copying of data but is also
        // mandatory to avoid conflicting read/write accesses on LAMA arrays

        LAMA_LOG_INFO( logger, typeName() << ": self assign, skipped, matrix = " << other )

        return;
    }

    LAMA_LOG_INFO( logger, typeName() << ": assign " << other )

    // Nearly the same routine as MatrixStorage::assign but here we
    // do not need any temporary data for ia, ja, and values

    other.buildCSRData( mIa, mJa, mValues );

    // actualize my member variables (class CSRStorage)

    _MatrixStorage::_assign( other ); // copy sizes

    mNumValues = mJa.size();

    mDiagonalProperty = checkDiagonalProperty();

    buildRowIndexes();

    check( "assign" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::assignTranspose( const MatrixStorage<ValueType>& other )
{
    LAMA_LOG_INFO( logger, *this << ": (CSR) assign transpose " << other )

    _MatrixStorage::_assignTranspose( other );

    // pass LAMAArrays of this storage to build the values in it

    if ( &other == this )
    {
        LAMAArray<IndexType> tmpIA;
        LAMAArray<IndexType> tmpJA;
        LAMAArray<ValueType> tmpValues;

        other.buildCSCData( tmpIA, tmpJA, tmpValues );
        swap( tmpIA, tmpJA, tmpValues );
    }
    else
    {
        other.buildCSCData( mIa, mJa, mValues );
        mNumValues = mJa.size();
        mDiagonalProperty = checkDiagonalProperty();
        buildRowIndexes();
    }

    // actualize my member variables (class CSRStorage)

    check( "assignTranspose" );

}

/* --------------------------------------------------------------------------- */

template<typename T>
void CSRStorage<T>::copyTo( _MatrixStorage& other ) const
{
    // Compressed sparse column data can be set directly to other matrix

    other.setCSRData( mNumRows, mNumColumns, mNumValues, mIa, mJa, mValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::buildCSR(
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr loc ) const
{
    ReadAccess<IndexType> inIA( mIa, loc );

    //build number of values per row into ia
    if ( ja == NULL || values == NULL )
    {
        WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows );

        LAMA_INTERFACE_FN( offsets2sizes, loc, CSRUtils, Offsets )
        LAMA_CONTEXT_ACCESS( loc )

        offsets2sizes( csrIA.get(), inIA.get(), mNumRows );
        return;
    }

    // copy the offset array ia and ja
    {
        ReadAccess<IndexType> inJA( mJa, loc );
        WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );
        WriteOnlyAccess<IndexType> csrJA( *ja, loc, mNumValues );

        LAMA_CONTEXT_ACCESS( loc )
        LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, IndexType, IndexType )
        set( csrIA.get(), inIA.get(), mNumRows + 1 );
        set( csrJA.get(), inJA.get(), mNumValues );
    }

    // copy values
    {
        ReadAccess<ValueType> inValues( mValues, loc );
        WriteOnlyAccess<OtherValueType> csrValues( *values, loc, mNumValues );

        LAMA_CONTEXT_ACCESS( loc )
        LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, OtherValueType, ValueType )
        set( csrValues.get(), inValues.get(), mNumValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildCSCData(
    LAMAArray<IndexType>& colIA,
    LAMAArray<IndexType>& colJA,
    LAMAArray<ValueType>& colValues ) const
{
    LAMA_LOG_INFO( logger, *this << ": buildCSCData by call of CSR2CSC" )

    // build the CSC data directly on the device where this matrix is located.

    this->convertCSR2CSC( colIA, colJA, colValues, mNumColumns, mIa, mJa, mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::splitHalo(
    MatrixStorage<ValueType>& localData,
    MatrixStorage<ValueType>& haloData,
    Halo& halo,
    const Distribution& colDist,
    const Distribution* rowDist ) const
{
    LAMA_REGION( "Storage.splitHalo" )

    LAMA_LOG_INFO( logger, *this << ": split CSR according to column distribution " << colDist )

    LAMA_ASSERT_EQUAL( mNumColumns, colDist.getGlobalSize() )

    if ( colDist.isReplicated() )
    {
        // if there is no column distribution, halo is not needed

        if ( rowDist )
        {
            localData.localize( *this, *rowDist );
        }
        else
        {
            localData.assign( *this );
        }

        haloData.allocate( mNumRows, 0 );

        halo = Halo(); // empty halo schedule

        return;
    }

    IndexType numRows = mNumRows;

    // check optional row distribution if specified

    if ( rowDist )
    {
        LAMA_LOG_INFO( logger, *this << ": split also localizes for " << *rowDist )
        LAMA_ASSERT_EQUAL( mNumRows, rowDist->getGlobalSize() )
        numRows = rowDist->getLocalSize();
    }

    LAMAArray<IndexType> localIA;
    LAMAArray<IndexType> localJA;
    LAMAArray<ValueType> localValues;

    LAMAArray<IndexType> haloIA;
    LAMAArray<IndexType> haloJA;
    LAMAArray<ValueType> haloValues;

    StorageMethods<ValueType>::splitCSR( localIA, localJA, localValues, haloIA, haloJA, haloValues, 
                                         mIa, mJa, mValues, colDist, rowDist );

    LAMA_ASSERT_EQUAL_DEBUG( localIA.size(), numRows + 1 )
    LAMA_ASSERT_EQUAL_DEBUG( haloIA.size(), numRows + 1 )

    const IndexType haloNumValues = haloJA.size();
    const IndexType localNumValues = localJA.size();

    LAMA_LOG_INFO( logger, *this << ": split into " << localNumValues << " local non-zeros "
                   " and " << haloNumValues << " halo non-zeros" )

    const IndexType localNumColumns = colDist.getLocalSize();

    IndexType haloNumColumns; // will be available after remap

    // build the halo by the non-local indexes

    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );

    LAMA_LOG_INFO( logger, "build halo: " << halo )

    localData.setCSRData( numRows, localNumColumns, localNumValues, localIA, localJA, localValues );

    localData.check( "local part after split" );

    // halo data is expected to have many empty rows, so enable compressing with row indexes

    haloData.setCompressThreshold( 0.5 );

    haloData.setCSRData( numRows, haloNumColumns, haloNumValues, haloIA, haloJA, haloValues );

    haloData.check( "halo part after split" );

    LAMA_LOG_INFO( logger,
                   "Result of split: local storage = " << localData << ", halo storage = " << haloData << ", halo = " << halo )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesVector(
    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    LAMA_REGION( "Storage.CSR.timesVector" )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( result.size(), mNumRows )

    if ( ( beta != 0.0 ) && ( result != y ) )
    {
        LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )
    }

    ContextPtr loc = getContextPtr();

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_T( sparseGEMV, loc, CSRUtils, Mult, ValueType )
    LAMA_INTERFACE_FN_T( normalGEMV, loc, CSRUtils, Mult, ValueType )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        if ( mRowIndexes.size() > 0 && ( beta == 1.0 ) )
        {
            // y += alpha * thisMatrix * x, can take advantage of row indexes

            IndexType numNonZeroRows = mRowIndexes.size();
            ReadAccess<IndexType> rows( mRowIndexes, loc );

            LAMA_CONTEXT_ACCESS( loc )
            sparseGEMV( wResult.get(), alpha, rX.get(), numNonZeroRows, rows.get(), csrIA.get(), csrJA.get(),
                        csrValues.get(), NULL );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            LAMA_CONTEXT_ACCESS( loc )
            normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumValues,
                        csrIA.get(), csrJA.get(), csrValues.get(), NULL );
        }
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc )
        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumValues,
                    csrIA.get(), csrJA.get(), csrValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::vectorTimesMatrix(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    LAMA_REGION( "Storage.CSR.VectorTimesMatrix" )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    LAMA_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if ( ( beta != 0.0 ) && ( result != y ) )
    {
        LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    ContextPtr loc = getContextPtr();

    LAMA_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )

    LAMA_INTERFACE_FN_T( sparseGEVM, loc, CSRUtils, Mult, ValueType )
    LAMA_INTERFACE_FN_T( normalGEVM, loc, CSRUtils, Mult, ValueType )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        if ( mRowIndexes.size() > 0 && ( beta == 1.0 ) )
        {
            // y += alpha * thisMatrix * x, can take advantage of row indexes

            IndexType numNonZeroRows = mRowIndexes.size();
            ReadAccess<IndexType> rows( mRowIndexes, loc );

            LAMA_CONTEXT_ACCESS( loc )
            sparseGEVM( wResult.get(), alpha, rX.get(), mNumColumns, numNonZeroRows, rows.get(), csrIA.get(), csrJA.get(),
                        csrValues.get(), NULL );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            LAMA_CONTEXT_ACCESS( loc )
            normalGEVM( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns,
                        csrIA.get(), csrJA.get(), csrValues.get(), NULL );
        }
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc )
        normalGEVM( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns,
                    csrIA.get(), csrJA.get(), csrValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesVectorN(
    LAMAArrayView<ValueType> result,
    const IndexType n,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    LAMA_REGION( "Storage.CSR.timesVectorN" )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), n * mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( result.size(), n * mNumRows )

    if ( ( beta != 0.0 ) && ( result != y ) )
    {
        LAMA_ASSERT_EQUAL_ERROR( y.size(), n * mNumRows )
    }

    ContextPtr loc = getContextPtr();

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVectorN on " << *loc )

    LAMA_INTERFACE_FN_T( gemm, loc, CSRUtils, Mult, ValueType )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        LAMA_CONTEXT_ACCESS( loc )
        gemm( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, n, mNumColumns, csrIA.get(), csrJA.get(),
              csrValues.get(), NULL );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc )
        gemm( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, n, mNumColumns, csrIA.get(), csrJA.get(),
              csrValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::matrixTimesVectorAsync(
    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": matrixTimesVectorAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    LAMA_REGION( "Storage.CSR.timesVectorAsync" )

    ContextPtr loc = getContextPtr();

    // Note: checks will be done by asynchronous task in any case
    //       and exception in tasks are handled correctly

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVectorAsync on " << *loc )

    if ( loc->getType() == Context::Host )
    {
        // execution as separate thread

        void (CSRStorage::*pf)(
            LAMAArrayView<ValueType>,
            const ValueType,
            const LAMAArrayConstView<ValueType>,
            const ValueType,
            const LAMAArrayConstView<ValueType> ) const

        = &CSRStorage<ValueType>::matrixTimesVector;

        using boost::bind;

        LAMA_LOG_INFO( logger, *this << ": matrixTimesVectorAsync on Host by own thread" )

		return new TaskSyncToken( bind( pf, this, result, alpha, x, beta, y ) );
    }

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( result.size(), mNumRows )

    if ( ( beta != 0.0 ) && ( result != y ) )
    {
        LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )
    }

    LAMA_INTERFACE_FN_T( sparseGEMV, loc, CSRUtils, Mult, ValueType )
    LAMA_INTERFACE_FN_T( normalGEMV, loc, CSRUtils, Mult, ValueType )

    auto_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > csrIA( new ReadAccess<IndexType>( mIa, loc ) );
    shared_ptr<ReadAccess<IndexType> > csrJA( new ReadAccess<IndexType>( mJa, loc ) );
    shared_ptr<ReadAccess<ValueType> > csrValues( new ReadAccess<ValueType>( mValues, loc ) );
    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        if ( mRowIndexes.size() > 0 && ( beta == 1.0 ) )
        {
            // y += alpha * thisMatrix * x, can take advantage of row indexes

            IndexType numNonZeroRows = mRowIndexes.size();

            shared_ptr<ReadAccess<IndexType> > rows( new ReadAccess<IndexType>( mRowIndexes, loc ) );

            syncToken->pushAccess( rows );

            LAMA_CONTEXT_ACCESS( loc )

            sparseGEMV( wResult->get(), alpha, rX->get(), numNonZeroRows, rows->get(), csrIA->get(), csrJA->get(),
                        csrValues->get(), syncToken.get() );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            LAMA_CONTEXT_ACCESS( loc )

            normalGEMV( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumColumns, mNumValues,
                        csrIA->get(), csrJA->get(), csrValues->get(), syncToken.get() );
        }

        syncToken->pushAccess( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumRows ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        LAMA_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumColumns, mNumValues,
                    csrIA->get(), csrJA->get(), csrValues->get(), syncToken.get() );

        syncToken->pushAccess( wResult );
        syncToken->pushAccess( rY );
    }

    syncToken->pushAccess( csrIA );
    syncToken->pushAccess( csrJA );
    syncToken->pushAccess( csrValues );
    syncToken->pushAccess( rX );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::vectorTimesMatrixAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    LAMA_REGION( "Storage.CSR.vectorTimesMatrixAsync" )

    ContextPtr loc = getContextPtr();

    // Note: checks will be done by asynchronous task in any case
    //       and exception in tasks are handled correctly

    LAMA_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on " << *loc )

    if ( loc->getType() == Context::Host )
    {
        // execution as separate thread

        void (CSRStorage::*pf)(
            LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>& ) const

        = &CSRStorage<ValueType>::vectorTimesMatrix;

        using boost::bind;

        LAMA_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on Host by own thread" )

        using boost::ref;

        return new TaskSyncToken( bind( pf, this, ref( result ), alpha, ref( x ), beta, ref( y ) ) );
    }

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    LAMA_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if ( ( beta != 0.0 ) && ( result != y ) )
    {
        LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    LAMA_INTERFACE_FN_T( sparseGEVM, loc, CSRUtils, Mult, ValueType )
    LAMA_INTERFACE_FN_T( normalGEVM, loc, CSRUtils, Mult, ValueType )

    auto_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > csrIA( new ReadAccess<IndexType>( mIa, loc ) );
    shared_ptr<ReadAccess<IndexType> > csrJA( new ReadAccess<IndexType>( mJa, loc ) );
    shared_ptr<ReadAccess<ValueType> > csrValues( new ReadAccess<ValueType>( mValues, loc ) );
    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        // only write access for y, no read access for result

        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        if ( mRowIndexes.size() > 0 && ( beta == 1.0 ) )
        {
            // y += alpha * thisMatrix * x, can take advantage of row indexes

            IndexType numNonZeroRows = mRowIndexes.size();

            shared_ptr<ReadAccess<IndexType> > rows( new ReadAccess<IndexType>( mRowIndexes, loc ) );

            syncToken->pushAccess( rows );

            LAMA_CONTEXT_ACCESS( loc )

            sparseGEVM( wResult->get(), alpha, rX->get(), mNumColumns, numNonZeroRows, rows->get(), csrIA->get(), csrJA->get(),
                        csrValues->get(), syncToken.get() );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            LAMA_CONTEXT_ACCESS( loc )

            normalGEVM( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumColumns,
                        csrIA->get(), csrJA->get(), csrValues->get(), syncToken.get() );
        }

        syncToken->pushAccess( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumColumns ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        LAMA_CONTEXT_ACCESS( loc )

        normalGEVM( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumColumns,
                    csrIA->get(), csrJA->get(), csrValues->get(), syncToken.get() );

        syncToken->pushAccess( wResult );
        syncToken->pushAccess( rY );
    }

    syncToken->pushAccess( csrIA );
    syncToken->pushAccess( csrJA );
    syncToken->pushAccess( csrValues );
    syncToken->pushAccess( rX );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.CSR.jacobiIterate" )

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

    LAMA_INTERFACE_FN_T( jacobi, loc, CSRUtils, Solver, ValueType )

    WriteAccess<ValueType> wSolution( solution, loc );
    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    // Due to diagonal property there is no advantage by taking row indexes

    LAMA_CONTEXT_ACCESS( loc )

    jacobi( wSolution.get(), csrIA.get(), csrJA.get(), csrValues.get(), rOldSolution.get(), rRhs.get(), omega, mNumRows,
            NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.CSR.jacobiIterateHalo" )

    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumRows() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumColumns() )
    LAMA_ASSERT_DEBUG( localStorage.hasDiagonalProperty(), localStorage << ": has not diagonal property" )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    const CSRStorage<ValueType>* csrLocal;

    if ( localStorage.getFormat() == CSR )
    {
        csrLocal = dynamic_cast<const CSRStorage<ValueType>*>( &localStorage );
        LAMA_ASSERT_DEBUG( csrLocal, "could not cast to CSRStorage " << localStorage )
    }
    else
    {
        // either copy localStorage to CSR (not recommended) or
        // just get the diagonal in localValues and set order in localIA
        LAMA_THROWEXCEPTION( "local stroage is not CSR" )
    }

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobiHalo, loc, CSRUtils, Solver, ValueType )

    {
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<IndexType> localIA( csrLocal->mIa, loc );
        ReadAccess<ValueType> localValues( csrLocal->mValues, loc );
        ReadAccess<IndexType> haloIA( mIa, loc );
        ReadAccess<IndexType> haloJA( mJa, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

        const IndexType numNonEmptyRows = mRowIndexes.size();

        LAMA_LOG_INFO( logger, "#row indexes = " << numNonEmptyRows )

        if ( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );

            LAMA_CONTEXT_ACCESS( loc )

            jacobiHalo( wSolution.get(), localIA.get(), localValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                        haloRowIndexes.get(), rOldHaloSolution.get(), omega, numNonEmptyRows );
        }
        else
        {
            LAMA_CONTEXT_ACCESS( loc )

            jacobiHalo( wSolution.get(), localIA.get(), localValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                        NULL, rOldHaloSolution.get(), omega, mNumRows );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const LAMAArray<ValueType>* localDiagonal,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.CSR.jacobiIterateHalo" )

    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobiHaloWithDiag, loc, CSRUtils, Solver, ValueType )

    {
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<ValueType> localDiagValues( *localDiagonal, loc );
        ReadAccess<IndexType> haloIA( mIa, loc );
        ReadAccess<IndexType> haloJA( mJa, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

        const IndexType numNonEmptyRows = mRowIndexes.size();

        LAMA_LOG_INFO( logger, "#row indexes = " << numNonEmptyRows )

        if ( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );

            LAMA_CONTEXT_ACCESS( loc )

            jacobiHaloWithDiag( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
				haloRowIndexes.get(), rOldHaloSolution.get(), omega, numNonEmptyRows );
        }
        else
        {
            LAMA_CONTEXT_ACCESS( loc )

            jacobiHaloWithDiag( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
				NULL, rOldHaloSolution.get(), omega, mNumRows );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const ValueType beta,
    const MatrixStorage<ValueType>& b )
{
    LAMA_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", with A = " << a << ", B = " << b )

    LAMA_REGION( "Storage.CSR.plusMatrix" )

    // a and b have to be CSR storages, otherwise create temporaries.

    const CSRStorage<ValueType>* csrA = NULL;
    const CSRStorage<ValueType>* csrB = NULL;

    // Define shared pointers in case we need temporaries

    boost::shared_ptr<CSRStorage<ValueType> > tmpA;
    boost::shared_ptr<CSRStorage<ValueType> > tmpB;

    if ( a.getFormat() == CSR )
    {
        csrA = dynamic_cast<const CSRStorage<ValueType>*>( &a );
        LAMA_ASSERT_DEBUG( csrA, "could not cast to CSRStorage " << a )
    }
    else
    {
        LAMA_UNSUPPORTED( a << ": will be converted to CSR for matrix multiply" )
        tmpA = boost::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( a ) );
        csrA = tmpA.get();
    }

    if ( b.getFormat() == CSR )
    {
        csrB = dynamic_cast<const CSRStorage<ValueType>*>( &b );
        LAMA_ASSERT_DEBUG( csrB, "could not cast to CSRStorage " << b )
    }
    else
    {
        LAMA_UNSUPPORTED( b << ": will be converted to CSR for matrix multiply" )
        tmpB = boost::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( b ) );
        csrB = tmpB.get();
    }

    // compute where target data will be

    ContextPtr loc = this->getContextPtr();

    matrixAddMatrixCSR( alpha, *csrA, beta, *csrB, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& c )
{
    LAMA_LOG_INFO( logger,
                   "this = " << alpha << " +* A * B + " << beta << " * C, with " << "A = " << a << ", B = " << b << ", C = " << c )

    LAMA_REGION( "Storage.CSR.timesMatrix" )

    // a and b have to be CSR storages, otherwise create temporaries.

    const CSRStorage<ValueType>* csrA = NULL;
    const CSRStorage<ValueType>* csrB = NULL;
    const CSRStorage<ValueType>* csrC = NULL;

    // Define two shared pointers in case we need temporaries

    boost::shared_ptr<CSRStorage<ValueType> > tmpA;
    boost::shared_ptr<CSRStorage<ValueType> > tmpB;
    boost::shared_ptr<CSRStorage<ValueType> > tmpC;

    if ( a.getFormat() == CSR )
    {
        csrA = dynamic_cast<const CSRStorage<ValueType>*>( &a );
        LAMA_ASSERT_DEBUG( csrA, "could not cast to CSRStorage " << a )
    }
    else
    {
        LAMA_UNSUPPORTED( a << ": will be converted to CSR for matrix multiply" )
        tmpA = boost::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( a ) );
        csrA = tmpA.get();
    }

    if ( b.getFormat() == CSR )
    {
        csrB = dynamic_cast<const CSRStorage<ValueType>*>( &b );
        LAMA_ASSERT_DEBUG( csrB, "could not cast to CSRStorage " << b )
    }
    else
    {
        LAMA_UNSUPPORTED( b << ": will be converted to CSR for matrix multiply" )
        tmpB = boost::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( b ) );
        csrB = tmpB.get();
    }

    if ( beta != 0.0 )
    {
        // c temporary needed if not correct format/type or aliased to this

        if ( ( c.getFormat() == CSR ) && ( &c != this ) )
        {
            csrC = dynamic_cast<const CSRStorage<ValueType>*>( &c );
            LAMA_ASSERT_DEBUG( csrC, "could not cast to CSRStorage " << c )
        }
        else
        {
            LAMA_UNSUPPORTED( c << ": CSR temporary required for matrix add" )
            tmpC = boost::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( c ) );
            csrC = tmpC.get();
        }

    }

    // now we have in any case all arguments as CSR Storage

    ContextPtr loc = ContextFactory::getContext( Context::Host );
    if ( a.getContext().getType() == b.getContext().getType() )
    {
        loc = a.getContextPtr();
    }

    ContextPtr saveContext = getContextPtr();

    CSRStorage<ValueType> tmp1;
    tmp1.matrixTimesMatrixCSR( alpha, *csrA, *csrB, loc );
    tmp1.setContext( loc );

    if ( beta != 0 )
    {
        CSRStorage<ValueType> tmp2;
        tmp2.matrixAddMatrixCSR( 1.0, tmp1, beta, *csrC, loc );
        swap( tmp2 );
    }
    else
    {
        swap( tmp1 );
    }

    this->setContext( saveContext );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixAddMatrixCSR(
    const ValueType alpha,
    const CSRStorage<ValueType>& a,
    const ValueType beta,
    const CSRStorage<ValueType>& b,
    const ContextPtr loc )
{
    LAMA_LOG_INFO( logger,
                   "this = " << alpha << " * A + " << beta << " * B, with " << "A = " << a << ", B = " << b << ", all are CSR" )

//    // TODO: just temporary, MAKE loc const again!
//    loc = ContextFactory::getContext( Context::Host );

    LAMA_INTERFACE_FN( matrixAddSizes, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_T( matrixAdd, loc, CSRUtils, Mult, ValueType )

    if ( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary
 
        CSRStorage<ValueType> tmp;

        tmp.matrixAddMatrixCSR( alpha, a, beta, b, loc );

        swap( tmp );  // safe as tmp will be destroyed afterwards

        return;
    }

    LAMA_REGION( "Storage.CSR.addMatrixCSR" )

    allocate( a.getNumRows(), a.getNumColumns() );

    LAMA_ASSERT_EQUAL_ERROR( mNumRows, b.getNumRows() )
    LAMA_ASSERT_EQUAL_ERROR( mNumColumns, b.getNumColumns() )

    mDiagonalProperty = ( mNumRows == mNumColumns );

    {
        ReadAccess<IndexType> aIa( a.getIA(), loc );
        ReadAccess<IndexType> aJa( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );

        ReadAccess<IndexType> bIa( b.getIA(), loc );
        ReadAccess<IndexType> bJa( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );

        // Step 1: compute row sizes of C, build offsets
        LAMA_LOG_DEBUG( logger, "Determing sizes of result matrix C" )

        WriteOnlyAccess<IndexType> cIa( mIa, loc, mNumRows + 1 );

        LAMA_CONTEXT_ACCESS( loc )

        mNumValues = matrixAddSizes( cIa.get(), mNumRows, mNumColumns, mDiagonalProperty, aIa.get(), aJa.get(),
                                     bIa.get(), bJa.get() );

        // Step 2: fill in ja, values

        LAMA_LOG_DEBUG( logger, "Compute the sparse values, # = " << mNumValues )

        WriteOnlyAccess<IndexType> cJa( mJa, loc, mNumValues );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValues );

        matrixAdd( cJa.get(), cValues.get(), cIa.get(), mNumRows, mNumColumns, mDiagonalProperty, alpha, aIa.get(),
                   aJa.get(), aValues.get(), beta, bIa.get(), bJa.get(), bValues.get() );
    }

    LAMA_LOG_DEBUG( logger, *this << ": compress by removing zero elements" )

    // Computation of C might have produced some zero elements

    //compress();

    check( "result of matrix + matrix" ); // just verify for a correct matrix
}

/* --------------------------------------------------------------------------- */

//TODO: just temporary
template<typename ValueType>
void CSRStorage<ValueType>::setNumValues( const IndexType numValues )
{
    mNumValues = numValues;
}

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesMatrixCSR(
    const ValueType alpha,
    const CSRStorage<ValueType>& a,
    const CSRStorage<ValueType>& b,
    const ContextPtr loc )
{
    LAMA_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B, with " << "A = " << a << ", B = " << b << ", all are CSR" << ", Context = " << loc->getType() )

    LAMA_INTERFACE_FN( matrixMultiplySizes, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_T( matrixMultiply, loc, CSRUtils, Mult, ValueType )

    LAMA_ASSERT_ERROR( &a != this, "matrixTimesMatrix: alias of a with this result matrix" )
    LAMA_ASSERT_ERROR( &b != this, "matrixTimesMatrix: alias of b with this result matrix" )

    LAMA_ASSERT_EQUAL_ERROR( a.getNumColumns(), b.getNumRows() )

    IndexType k = a.getNumColumns();

    LAMA_REGION( "Storage.CSR.timesMatrixCSR" )

    allocate( a.getNumRows(), b.getNumColumns() );

    mDiagonalProperty = ( mNumRows == mNumColumns );

    {
        ReadAccess<IndexType> aIA( a.getIA(), loc );
        ReadAccess<IndexType> aJA( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );

        ReadAccess<IndexType> bIA( b.getIA(), loc );
        ReadAccess<IndexType> bJA( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );

        WriteOnlyAccess<IndexType> cIA( mIa, loc, mNumRows + 1 );

        LAMA_CONTEXT_ACCESS( loc )

        mNumValues = matrixMultiplySizes( cIA.get(), mNumRows, mNumColumns, k, mDiagonalProperty, aIA.get(), aJA.get(),
                                          bIA.get(), bJA.get() );

        WriteOnlyAccess<IndexType> cJa( mJa, loc, mNumValues );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValues );

        matrixMultiply( cIA.get(), cJa.get(), cValues.get(), 
                        mNumRows, mNumColumns, k, alpha, mDiagonalProperty, aIA.get(),
                        aJA.get(), aValues.get(), bIA.get(), bJA.get(), bValues.get() );
    }

    // TODO: check this!
//    compress();
    buildRowIndexes();
//    check( "result of matrix x matrix" ); // just verify for a correct matrix
//    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::maxNorm() const
{
    // no more checks needed here

    LAMA_LOG_INFO( logger, *this << ": maxNorm()" )

    if ( mNumValues == 0 )
    {
        return 0.0f;
    }

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, Utils, Reductions, ValueType )

    ReadAccess<ValueType> csrValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal( csrValues.get(), mNumValues );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    LAMA_REGION( "Storage.CSR.maxDiffNorm" )

    LAMA_ASSERT_EQUAL_ERROR( mNumRows, other.getNumRows() )
    LAMA_ASSERT_EQUAL_ERROR( mNumColumns, other.getNumColumns() )

    LAMA_LOG_INFO( logger, *this << ": maxDiffNorm( " << other << " )" )

    boost::shared_ptr<CSRStorage<ValueType> > tmpOtherCSR;

    const CSRStorage<ValueType>* otherCSR;

    if ( other.getValueType() == this->getValueType() && ( other.getFormat() == CSR ) )
    {
        otherCSR = dynamic_cast<const CSRStorage<ValueType>*>( &other );
        LAMA_ASSERT_ERROR( otherCSR, other << ": could not cast to " << typeName() )
    }
    else
    {
        LAMA_UNSUPPORTED( other << ": converted to " << typeName() << " for maxDiffNorm" )
        tmpOtherCSR.reset( new CSRStorage<ValueType>( other ) );
        otherCSR = tmpOtherCSR.get();
    }

    return maxDiffNormImpl( *otherCSR );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::maxDiffNormImpl( const CSRStorage<ValueType>& other ) const
{
    // no more checks needed here

    LAMA_LOG_INFO( logger, *this << ": maxDiffNormImpl( " << other << " )" )

    if ( mNumRows == 0 )
    {
        return 0.0f;
    }

    ContextPtr loc = getContextPtr();

    bool sorted = mSortedRows && other.mSortedRows && ( mDiagonalProperty == other.mDiagonalProperty );

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxDiffVal, loc, CSRUtils, Reductions, ValueType )

    ReadAccess<IndexType> csrIA1( mIa, loc );
    ReadAccess<IndexType> csrJA1( mJa, loc );
    ReadAccess<ValueType> csrValues1( mValues, loc );

    ReadAccess<IndexType> csrIA2( other.mIa, loc );
    ReadAccess<IndexType> csrJA2( other.mJa, loc );
    ReadAccess<ValueType> csrValues2( other.mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxDiffVal( mNumRows, sorted, csrIA1.get(), csrJA1.get(), csrValues1.get(), csrIA2.get(),
                                      csrJA2.get(), csrValues2.get() );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>* CSRStorage<ValueType>::create() const
{
    return new CSRStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>* CSRStorage<ValueType>::copy() const
{
    return new CSRStorage<ValueType>( *this );
}

template<typename ValueType>
void CSRStorage<ValueType>::buildSparseRowSizes( LAMAArray<IndexType>& rowSizes ) const
{
    LAMA_LOG_DEBUG( logger, "copy nnz for each row in LAMAArray" );

    HostWriteOnlyAccess<IndexType> writeRowSizes( rowSizes, mNumRows );
    HostReadAccess<IndexType> csrIA( mIa );

    OpenMPCSRUtils::offsets2sizes( writeRowSizes.get(), csrIA.get(), mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildSparseRowData( LAMAArray<IndexType>& sparseJA,
                                                LAMAArray<ValueType>& sparseValues ) const
{
    LAMA_LOG_INFO( logger, *this << ": build sparse row data" );

    // for CSR format we can just copy arrays with column indexes and data values

    sparseJA     = mJa;
    sparseValues = mValues;
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_CSR_STORAGE_INSTANTIATE(z, I, _)                              \
template<>                                                                 \
const char* CSRStorage<ARITHMETIC_TYPE##I>::typeName()                     \
{                                                                          \
    return "CSRStorage<ARITHMETIC_TYPE##I>";                               \
}                                                                          \
                                                                           \
template class LAMA_DLL_IMPORTEXPORT CSRStorage<ARITHMETIC_TYPE##I> ;      \
                                                                           \
template void CSRStorage<ARITHMETIC_TYPE##I>::setCSRDataSwap(              \
    const IndexType numRows,                                               \
    const IndexType numColumns,                                            \
    const IndexType numValues,                                             \
    LAMAArray<IndexType>& ia,                                              \
    LAMAArray<IndexType>& ja,                                              \
    LAMAArray<ARITHMETIC_TYPE##I>& values,                                \
    const ContextPtr loc );                                                \

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_CSR_STORAGE_INSTANTIATE, _ )

#undef LAMA_CSR_STORAGE_INSTANTIATE

}
