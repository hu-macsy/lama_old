/**
 * @file CSRStorage.cpp
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
 * @brief Implementation and instantiation for template class CSRStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/CSRStorage.hpp>

// local library
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>


#include <scai/lama/storage/StorageMethods.hpp>

#include <scai/dmemo/Redistributor.hpp>


// internal scai libraries
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/assert.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/tasking/NoSyncToken.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;
using namespace dmemo;
using namespace utilskernel;

using common::unique_ptr;
using common::shared_ptr;
using common::TypeTraits;

using sparsekernel::CSRKernelTrait;
using sparsekernel::OpenMPCSRUtils;

using tasking::SyncToken;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, CSRStorage<ValueType>::logger, "MatrixStorage.CSRStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage()
    : CRTPMatrixStorage<CSRStorage<ValueType>,ValueType>( 0, 0 ), mNumValues( 0 ), mSortedRows( false )
{
    allocate( 0, 0 ); // creates at least mIa
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>::CSRStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )

    : CRTPMatrixStorage<CSRStorage<ValueType>,ValueType>()
{
    this->setCSRData( numRows, numColumns, numValues, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;

    stream << "CSRStorage " << mNumRows << " x " << mNumColumns << ", #values = " << mNumValues << endl;

    ContextPtr host = Context::getHostPtr();

    ReadAccess<IndexType> ia( mIa, host );
    ReadAccess<IndexType> ja( mJa, host );
    ReadAccess<ValueType> values( mValues, host );

    for( IndexType i = 0; i < mNumRows; i++ )
    {
        stream << "Row " << i << " ( " << ia[i] << " - " << ia[i + 1] << " ) :";

        for( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
        {
            stream << " " << ja[jj] << ":" << values[jj];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::clear()
{
    allocate( 0, 0 );   // sets everything correctly, ia array has one entry
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format::MatrixStorageFormat CSRStorage<ValueType>::getFormat() const
{
    return Format::CSR;
}

/* --------------------------------------------------------------------------- */

#ifndef SCAI_ASSERT_LEVEL_DEBUG
template<typename ValueType>
void CSRStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void CSRStorage<ValueType>::check( const char* msg ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows + 1, mIa.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mJa.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )

    // check ascending values in offset array mIa

    {
        static LAMAKernel<UtilKernelTrait::isSorted<IndexType> > isSorted;
        static LAMAKernel<UtilKernelTrait::getValue<IndexType> > getValue;

        ContextPtr loc = this->getContextPtr();
        isSorted.getSupportedContext( loc, getValue );

        ReadAccess<IndexType> csrIA( mIa, loc );

        SCAI_CONTEXT_ACCESS( loc )

        bool ascending = true; // check for ascending

        IndexType numValues = getValue[ loc ]( csrIA.get(), mNumRows );

        SCAI_ASSERT_ERROR(
            numValues == mNumValues,
            "ia[" << mNumRows << "] = " << numValues << ", expected " << mNumValues << ", msg = " << msg )

        SCAI_ASSERT_ERROR( isSorted[ loc ]( csrIA.get(), mNumRows + 1, ascending ),
                           *this << " @ " << msg << ": IA is illegal offset array" )
    }

    // check column indexes in JA

    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;

        ContextPtr loc = this->getContextPtr();
        validIndexes.getSupportedContext( loc );

        ReadAccess<IndexType> rJA( mJa, loc );

        SCAI_CONTEXT_ACCESS( loc )

        SCAI_ASSERT_ERROR( validIndexes[loc]( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool CSRStorage<ValueType>::checkDiagonalProperty() const
{
    // diagonal property is given if size of matrix is 0

    if( mNumRows == 0 || mNumColumns == 0 )
    {
        return true;
    }

    // non-zero sized matrix with no values has not diagonal property

    if( mNumValues == 0 )
    {
        return false;
    }

    static LAMAKernel<CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;

    ContextPtr loc = this->getContextPtr();
    hasDiagonalProperty.getSupportedContext( loc );

    //get read access
    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    SCAI_CONTEXT_ACCESS( loc )

    IndexType numDiagonals = std::min( mNumRows, mNumColumns );

    bool diagonalProperty = hasDiagonalProperty[loc]( numDiagonals, csrIA.get(), csrJA.get() );

    SCAI_LOG_DEBUG( logger, *this << ": diagonalProperty = " << diagonalProperty );

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_DEBUG( logger, "set identity, size = " << size )

    _MatrixStorage::setDimension( size, size );

    mNumValues = mNumRows;

    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    {
        ContextPtr loc = this->getContextPtr();
        setOrder.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> ia( mIa, loc, mNumRows + 1 );
        WriteOnlyAccess<IndexType> ja( mJa, loc, mNumValues );

        setOrder[ loc ]( ia.get(), mNumRows + 1 );
        setOrder[ loc ]( ja.get(), mNumRows );
    }

    {
        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<ValueType> values( mValues, loc, mNumValues );
        setVal[loc]( values.get(), mNumRows, ValueType( 1 ), utilskernel::reduction::COPY );
    }


    mDiagonalProperty = true; // obviously given for identity matrix
    mSortedRows = true; // obviously given for identity matrix

    // Note: we do not build row indexes, no row is empty

    SCAI_LOG_INFO( logger, *this << ": identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr /* loc */)
{
    ContextPtr loc = this->getContextPtr();

    if( ia.size() == numRows )
    {
        IndexType sumIA = HArrayUtils::reduce( ia, utilskernel::reduction::ADD );

        SCAI_ASSERT_EQUAL( numValues, sumIA, "sizes do not sum up to numValues" );
    }
    else if ( ia.size() == numRows + 1 )
    {
        bool ascending = true; // check increasing, ia[i] <= ia[i+1]

        SCAI_ASSERT( HArrayUtils::isSorted( ia, ascending ), "ia is invalid offset array, entries not ascending" )

        SCAI_ASSERT_EQUAL( numValues, HArrayUtils::getValImpl( ia, numRows ), "last entry in offsets must be numValues" );
    }
    else
    {
        COMMON_THROWEXCEPTION( "ia array with size = " << ia.size() << " illegal, #rows = " << numRows )
    }

    SCAI_ASSERT_EQUAL_ERROR( numValues, ja.size() );
    SCAI_ASSERT_EQUAL_ERROR( numValues, values.size() );

    SCAI_ASSERT( HArrayUtils::validIndexes( ja, numColumns ), "invalid column indexes, #cols = " << numColumns );

    // now we can copy all data

    mNumRows = numRows;
    mNumColumns = numColumns;
    mNumValues = numValues;

    SCAI_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    // storage data will be directly allocated on the location

    if( ia.size() == numRows )
    {
        {
            // reserve enough memory for mIa

            WriteOnlyAccess<IndexType> myIA( mIa, loc, mNumRows + 1 );
        }

        HArrayUtils::assign( mIa, ia, loc );
        HArrayUtils::scan( mIa, loc );
    }
    else
    {
        HArrayUtils::assign( mIa, ia, loc );
    }

    HArrayUtils::assign( mValues, values, loc );
    HArrayUtils::assign( mJa, ja, loc );

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
        ReadAccess<IndexType> csrIA( mIa );
        WriteAccess<IndexType> csrJA( mJa );
        WriteAccess<ValueType> csrValues( mValues );

        OpenMPCSRUtils::sortRowElements( csrJA.get(), csrValues.get(), csrIA.get(), mNumRows, diagonalProperty );
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
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<OtherValueType>& values,
    const ContextPtr /* loc */)
{
    //set necessary information
    mNumRows = numRows;
    mNumColumns = numColumns;
    mNumValues = numValues;

    SCAI_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    //swap arrays

    mIa.swap( ia );
    mJa.swap( ja );

    if( common::TypeTraits<ValueType>::stype == common::TypeTraits<OtherValueType>::stype)
    {
        mValues.swap( reinterpret_cast<HArray<ValueType>&>( values ) );
    }
    else
    {
        COMMON_THROWEXCEPTION( "ValueType mismatch" )
    }

    mDiagonalProperty = checkDiagonalProperty();

    // this builds only row indices if context is on host

    buildRowIndexes();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildRowIndexes()
{
    mRowIndexes.clear();

    if ( mDiagonalProperty )
    {
        SCAI_LOG_INFO( logger, "buildRowIndexes not done as diagonal property implies at least one entry per row" );
        return;
    }

    if( mNumRows == 0 )
    {
        return;
    }

    if( getContextPtr()->getType() != Context::Host )
    {
        SCAI_LOG_INFO( logger, "CSRStorage: build row indices is currently only implemented on host" )
    }

    // This routine is only available on the Host

    ContextPtr loc = Context::getHostPtr();

    ReadAccess<IndexType> csrIA( mIa, loc );

    IndexType nonZeroRows = OpenMPCSRUtils::countNonEmptyRowsByOffsets( csrIA.get(), mNumRows );

    float usage = float( nonZeroRows ) / float( mNumRows );

    if( usage >= mCompressThreshold )
    {
        SCAI_LOG_INFO( logger, "CSRStorage: do not build row indexes, usage = " << usage
                       << ", threshold = " << mCompressThreshold )
        return;
    }

    SCAI_LOG_INFO( logger, "CSRStorage: build row indexes, #entries = " << nonZeroRows )

    WriteOnlyAccess<IndexType> rowIndexes( mRowIndexes, loc, nonZeroRows );

    OpenMPCSRUtils::setNonEmptyRowsByOffsets( rowIndexes.get(), nonZeroRows, csrIA.get(), mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::redistributeCSR( const CSRStorage<ValueType>& other, const dmemo::Redistributor& redistributor )
{
    SCAI_REGION( "Storage.redistributeCSR" )

    const dmemo::Distribution& sourceDistribution = *redistributor.getSourceDistributionPtr();
    const dmemo::Distribution& targetDistribution = *redistributor.getTargetDistributionPtr();

    SCAI_LOG_INFO( logger,
                   other << ": redistribute of CSR<" << other.getValueType() << "> to CSR<" << this->getValueType() << " via " << redistributor )

    bool sameDist = false;

    // check for same distribution, either equal or both replicated

    if( sourceDistribution.isReplicated() && targetDistribution.isReplicated() )
    {
        sameDist = true;
    }
    else if( &sourceDistribution == &targetDistribution )
    {
        sameDist = true;
    }

    if( sameDist )
    {
        SCAI_LOG_INFO( logger, "redistributor with same source/target distribution" )

        assign( other );

        return; // so we are done
    }

    // check that source distribution fits with storage

    SCAI_ASSERT_EQUAL_ERROR( other.getNumRows(), sourceDistribution.getLocalSize() )

    if( &other == this )
    {
        // due to alias we need temporary array

        HArray<IndexType> targetIA;
        HArray<IndexType> targetJA;
        HArray<ValueType> targetValues;

        StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues, other.getIA(), other.getJA(),
                other.getValues(), redistributor );

        // we can swap the new arrays

        mIa.swap( targetIA );
        mJa.swap( targetJA );
        mValues.swap( targetValues );
    }
    else
    {
        StorageMethods<ValueType>::redistributeCSR( mIa, mJa, mValues, other.getIA(), other.getJA(), other.getValues(),
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
    SCAI_LOG_DEBUG( logger,
                    "~CSRStorage, size = " << mNumRows << " x " << mNumColumns << ", # non-zeros = " << mNumValues )
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
    // delete all old values

    mIa.purge();
    mJa.purge();
    mValues.purge();

    allocate( 0, 0 );   // sets everything correctly, ia array has one entry
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger,
                   "allocate CSR sparse matrix of size " << numRows << " x " << numColumns << ", numValues = 0" )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValues = 0;

    mJa.clear();
    mValues.clear();

    WriteOnlyAccess<IndexType> ia( mIa, mNumRows + 1 );

    // make a correct initialization for the offset array

    OpenMPUtils::setVal( ia.get(), mNumRows + 1, IndexType( 0 ), utilskernel::reduction::COPY  );

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::compress( const ValueType eps )
{
    static LAMAKernel<CSRKernelTrait::countNonZeros<ValueType> > countNonZeros;

    LArray<IndexType> newIa;

    {
        ContextPtr loc = this->getContextPtr();
        countNonZeros.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> ia( mIa, loc );
        ReadAccess<IndexType> ja( mJa, loc );
        ReadAccess<ValueType> values( mValues, loc );
        WriteOnlyAccess<IndexType> new_ia( newIa, loc, mNumRows + 1 );

        countNonZeros[loc]( new_ia.get(), ia.get(), ja.get(), values.get(), mNumRows, eps, mDiagonalProperty );
    }

    // now compute the new offsets from the sizes, gives also new numValues

    IndexType newNumValues = HArrayUtils::scan( newIa, this->getContextPtr() );

    SCAI_LOG_INFO( logger, "compress: " << newNumValues << " non-diagonal zero elements" )

    // ready if there are no new non-zero values

    if ( newNumValues == mNumValues )
    {
        return;
    }

    // All information is available how to fill the compressed data

    LArray<ValueType> newValues;
    LArray<IndexType> newJa;

    {
        static LAMAKernel<CSRKernelTrait::compress<ValueType> > compressData;

        ContextPtr loc = this->getContextPtr();
        compressData.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> new_ia( newIa, loc );

        ReadAccess<IndexType> ia( mIa, loc );
        ReadAccess<IndexType> ja( mJa, loc );
        ReadAccess<ValueType> values( mValues, loc );

        WriteOnlyAccess<IndexType> new_ja( newJa, newNumValues );
        WriteOnlyAccess<ValueType> new_values( newValues, newNumValues );

        compressData[loc]( new_ja.get(), new_values.get(), new_ia.get(),
                           ia.get(), ja.get(), values.get(), mNumRows,
                           eps, mDiagonalProperty );
    }

    // now switch in place to the new data

    mIa.swap( newIa );
    mJa.swap( newJa );
    mValues.swap( newValues );
    mNumValues = newNumValues;

    // Note: temporary data is freed implicitly
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
void CSRStorage<ValueType>::swap( HArray<IndexType>& ia, HArray<IndexType>& ja, HArray<ValueType>& values )
{
    SCAI_ASSERT_EQUAL_ERROR( ia.size(), mNumRows + 1 )

    IndexType numValues = 0;

    {
        ReadAccess<IndexType> csrIA( ia );
        numValues = csrIA[mNumRows];
    }

    SCAI_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, values.size() )

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
    stream << "CSRStorage<" << common::getScalarType<ValueType>() << ">("
           << " size = " << mNumRows << " x " << mNumColumns
           << ", nnz = " << mNumValues << ", diag = " << mDiagonalProperty << ", sorted = " << mSortedRows << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )
    const ReadAccess<IndexType> ia( mIa );
    const ReadAccess<IndexType> ja( mJa );
    const ReadAccess<ValueType> values( mValues );
    ValueType myValue = static_cast<ValueType>(0.0);

    SCAI_LOG_TRACE( logger, "search column in ja from " << ia[i] << ":" << ia[i + 1] )

    for( IndexType jj = ia[i]; jj < ia[i + 1]; ++jj )
    {
        IndexType col = ja[jj];

        SCAI_ASSERT_DEBUG( 0 <= col && col < mNumColumns,
                           "column index at pos " << jj << " = " << col << " out of range" )

        if( col == j )
        {
            SCAI_LOG_TRACE( logger, "found column j = " << j << " at " << jj << ", value = " << values[jj] )
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
LArray<IndexType>& CSRStorage<ValueType>::getIA()
{
    return mIa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
LArray<IndexType>& CSRStorage<ValueType>::getJA()
{
    return mJa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
LArray<ValueType>& CSRStorage<ValueType>::getValues()
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LArray<IndexType>& CSRStorage<ValueType>::getIA() const
{
    return mIa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LArray<IndexType>& CSRStorage<ValueType>::getJA() const
{
    return mJa;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LArray<ValueType>& CSRStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    if( !mDiagonalProperty )
    {
        COMMON_THROWEXCEPTION( "setDiagonal: matrix storage has not diagonal property." )
    }

    ReadAccess<IndexType> wIa( mIa );
    WriteAccess<ValueType> wValues( mValues );

    for( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wValues[wIa[i]] = value;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::setDiagonalImpl( const HArray<OtherValueType>& diagonal )
{
    IndexType numDiagonalElements = diagonal.size();

    if ( numDiagonalElements > mNumRows )
    {
        numDiagonalElements = mNumRows;
    }

    {
        static LAMAKernel<UtilKernelTrait::setScatter<ValueType, OtherValueType> > setScatter;

        ContextPtr loc = this->getContextPtr();
        setScatter.getSupportedContext( loc );

        SCAI_LOG_INFO( logger, "set diagonal<" << TypeTraits<ValueType>::id() << ", "
                       << TypeTraits<OtherValueType>::id() << "> ( " << numDiagonalElements << " ) @ " << *loc )

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
        ReadAccess<IndexType> csrIA( mIa, loc );

        WriteAccess<ValueType> wValues( mValues, loc );     // partial setting

        //  wValues[ wIa[ i ] ] = rDiagonal[ i ];

        setScatter[loc]( wValues.get(), csrIA.get(), rDiagonal.get(), numDiagonalElements );
    }

    if ( SCAI_LOG_TRACE_ON( logger ) )
    {
        SCAI_LOG_TRACE( logger, "CSR after setDiagonal" )
        print( std::cout );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void CSRStorage<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    IndexType n1;     // offset for row i in ja, values
    IndexType nrow;   // number of nonzero entries in row i

    {
        static LAMAKernel<UtilKernelTrait::getValue<IndexType> > getValue;

        // get ia[i], ia[i+1] from any location with valid values

        ContextPtr loc = mIa.getValidContext();
        getValue.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> ia( mIa, loc );

        n1 = getValue[loc]( ia.get(), i );
        nrow = getValue[loc]( ia.get(), i + 1 ) - n1;
    }

    static LAMAKernel<UtilKernelTrait::setVal<OtherType> > setVal;
    static LAMAKernel<UtilKernelTrait::setScatter<OtherType, ValueType> > setScatter;

    /// ContextPtr loc = Context::getHostPtr();

    ContextPtr loc = this->getContextPtr();
    setVal.getSupportedContext( loc, setScatter );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<OtherType> wRow( row, loc, mNumColumns );

    setVal[loc]    ( wRow.get(), mNumColumns, OtherType( 0 ), utilskernel::reduction::COPY );

    const ReadAccess<IndexType> ja( mJa, loc );
    const ReadAccess<ValueType> values( mValues, loc );

    setScatter[loc]( wRow.get(), ja.get() + n1, values.get() + n1, nrow );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::getDiagonalImpl( HArray<OtherValueType>& diagonal ) const
{
    const IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static LAMAKernel<UtilKernelTrait::setGather<OtherValueType, ValueType> > setGather;

    ContextPtr loc = this->getContextPtr();
    setGather.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<OtherValueType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<ValueType> rValues( mValues, loc );

    setGather[loc]( wDiagonal.get(), rValues.get(), csrIA.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::scaleImpl( const ValueType value )
{
    SCAI_ASSERT_EQUAL( mValues.size(), mNumValues, "size mismatch" )

    HArrayUtils::scale( mValues, value, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::conj()
{
    HArrayUtils::conj( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::scaleImpl( const HArray<OtherValueType>& diagonal )
{
    IndexType n = std::min( mNumRows, diagonal.size() );

    static LAMAKernel<CSRKernelTrait::scaleRows<ValueType, OtherValueType> > scaleRows;

    ContextPtr loc = this->getContextPtr();
    scaleRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    {
        ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
        ReadAccess<IndexType> csrIA( mIa, loc );
        WriteAccess<ValueType> csrValues( mValues, loc ); // updateAccess

        scaleRows[loc]( csrValues.get(), csrIA.get(), n, rDiagonal.get() );
    }

    if( SCAI_LOG_TRACE_ON( logger ) )
    {
        SCAI_LOG_TRACE( logger, "CSR after scale diagonal" )
        print( std::cout );
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
    if( &other == this )
    {
        // this special case avoids copying of data but is also
        // mandatory to avoid conflicting read/write accesses on LAMA arrays

        SCAI_LOG_INFO( logger, typeName() << ": self assign, skipped, matrix = " << other )

        return;
    }

    SCAI_LOG_INFO( logger, typeName() << ": assign " << other )

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
    SCAI_LOG_INFO( logger, *this << ": (CSR) assign transpose " << other )

    // pass HArrays of this storage to build the values in it

    if ( &other == this )
    {
        SCAI_LOG_INFO( logger, *this << ": (CSR) assign transpose in place" )

        HArray<IndexType> tmpIA;
        HArray<IndexType> tmpJA;
        HArray<ValueType> tmpValues;

        // do not change sizes before building CSC data

        other.buildCSCData( tmpIA, tmpJA, tmpValues );

        // sizes must be set correctly BEFORE swap

        _MatrixStorage::_assignTranspose( other );

        swap( tmpIA, tmpJA, tmpValues );  // sets all other data correctly
    }
    else
    {
        _MatrixStorage::_assignTranspose( other );

        SCAI_LOG_INFO( logger, *this << ": (CSR) assign transpose " << other )

        other.buildCSCData( mIa, mJa, mValues );

        mNumValues = mJa.size();
        mDiagonalProperty = checkDiagonalProperty();
        buildRowIndexes();
    }

    // actualize my member variables (class CSRStorage)

    check( "assignTranspose" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::copyTo( _MatrixStorage& other ) const
{
    // Compressed sparse column data can be set directly to other matrix

    other.setCSRData( mNumRows, mNumColumns, mNumValues, mIa, mJa, mValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>* CSRStorage<ValueType>::newMatrixStorage() const
{
    common::unique_ptr<CSRStorage<ValueType> > storage( new CSRStorage<ValueType>() );
    storage->setContextPtr( this->getContextPtr() );
    return storage.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void CSRStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<OtherValueType>* values,
    const ContextPtr prefLoc ) const
{
    // in case of ja == NULL we need only the size array, build it from offsets

    if ( ja == NULL || values == NULL )
    {
        static LAMAKernel<CSRKernelTrait::offsets2sizes> offsets2sizes;

        ContextPtr loc = prefLoc;
        offsets2sizes.getSupportedContext( loc );

        ReadAccess<IndexType> inIA( mIa, loc );
        WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        offsets2sizes[ loc ]( csrIA.get(), inIA.get(), mNumRows );

        return;
    }

    // copy the offset array ia and ja

    {
        static LAMAKernel<UtilKernelTrait::set<IndexType, IndexType> > setIndexes;

        ContextPtr loc = prefLoc;
        setIndexes.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> inIA( mIa, loc );
        ReadAccess<IndexType> inJA( mJa, loc );
        WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );
        WriteOnlyAccess<IndexType> csrJA( *ja, loc, mNumValues );

        setIndexes[ loc ]( csrIA.get(), inIA.get(), mNumRows + 1, utilskernel::reduction::COPY );
        setIndexes[ loc ]( csrJA.get(), inJA.get(), mNumValues, utilskernel::reduction::COPY );
    }

    // copy values
    {
        static LAMAKernel<UtilKernelTrait::set<OtherValueType, ValueType> > setValues;

        ContextPtr loc = prefLoc;
        setValues.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<ValueType> inValues( mValues, loc );
        WriteOnlyAccess<OtherValueType> csrValues( *values, loc, mNumValues );

        setValues[ loc ]( csrValues.get(), inValues.get(), mNumValues, utilskernel::reduction::COPY );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildCSCData(
    HArray<IndexType>& colIA,
    HArray<IndexType>& colJA,
    HArray<ValueType>& colValues ) const
{
    SCAI_LOG_INFO( logger, *this << ": buildCSCData by call of CSR2CSC" )

    // build the CSC data directly on the device where this matrix is located.

    this->convertCSR2CSC( colIA, colJA, colValues, mNumColumns, mIa, mJa, mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::splitHalo(
    MatrixStorage<ValueType>& localData,
    MatrixStorage<ValueType>& haloData,
    dmemo::Halo& halo,
    const dmemo::Distribution& colDist,
    const dmemo::Distribution* rowDist ) const
{
    SCAI_REGION( "Storage.splitHalo" )

    SCAI_LOG_INFO( logger, *this << ": split CSR according to column distribution " << colDist )

    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, colDist.getGlobalSize() )

    if( colDist.isReplicated() )
    {
        // if there is no column distribution, halo is not needed

        if( rowDist )
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

    if( rowDist )
    {
        SCAI_LOG_INFO( logger, *this << ": split also localizes for " << *rowDist )
        SCAI_ASSERT_EQUAL_ERROR( mNumRows, rowDist->getGlobalSize() )
        numRows = rowDist->getLocalSize();
    }

    HArray<IndexType> localIA;
    HArray<IndexType> localJA;
    HArray<ValueType> localValues;

    HArray<IndexType> haloIA;
    HArray<IndexType> haloJA;
    HArray<ValueType> haloValues;

    StorageMethods<ValueType>::splitCSR( localIA, localJA, localValues, haloIA, haloJA, haloValues, mIa, mJa, mValues,
                                         colDist, rowDist );

    SCAI_ASSERT_EQUAL_DEBUG( localIA.size(), numRows + 1 )
    SCAI_ASSERT_EQUAL_DEBUG( haloIA.size(), numRows + 1 )

    const IndexType haloNumValues = haloJA.size();
    const IndexType localNumValues = localJA.size();

    SCAI_LOG_INFO( logger,
                   *this << ": split into " << localNumValues << " local non-zeros " " and " << haloNumValues << " halo non-zeros" )

    const IndexType localNumColumns = colDist.getLocalSize();

    IndexType haloNumColumns; // will be available after remap

    // build the halo by the non-local indexes

    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );

    SCAI_LOG_INFO( logger, "build halo: " << halo )

    localData.setCSRData( numRows, localNumColumns, localNumValues, localIA, localJA, localValues );

    localData.check( "local part after split" );

    // halo data is expected to have many empty rows, so enable compressing with row indexes

    haloData.setCompressThreshold( 0.5 );

    haloData.setCSRData( numRows, haloNumColumns, haloNumValues, haloIA, haloJA, haloValues );

    haloData.check( "halo part after split" );

    SCAI_LOG_INFO( logger,
                   "Result of split: local storage = " << localData << ", halo storage = " << haloData << ", halo = " << halo )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    bool async = false; // synchronously execution, no SyncToken required

    SyncToken* token = gemv( result, alpha, x, beta, y, async );

    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    bool async = false;  // synchronously execution, no SyncToken required

    SyncToken* token = gevm( result, alpha, x, beta, y, async );

    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixTimesVectorN(
    HArray<ValueType>& result,
    const IndexType n,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.CSR.timesVectorN" )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), n * mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), n * mNumRows )

    if( ( beta != common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), n * mNumRows )
    }

    static LAMAKernel<CSRKernelTrait::gemm<ValueType> > gemm;

    ContextPtr loc = this->getContextPtr();
    gemm.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVectorN on " << *loc )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    // due to possible alias of result and y, write access must follow read(y)

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

    SCAI_CONTEXT_ACCESS( loc )

    gemm[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, n, mNumColumns,
               csrIA.get(), csrJA.get(), csrValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::sparseGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    static LAMAKernel<CSRKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = this->getContextPtr();
    sparseGEMV.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteAccess<ValueType> wResult( result, loc );

    // result += alpha * thisMatrix * x, can take advantage of row indexes

    IndexType numNonZeroRows = mRowIndexes.size();

    ReadAccess<IndexType> rows( mRowIndexes, loc );

    sparseGEMV[loc]( wResult.get(), alpha, rX.get(), numNonZeroRows, rows.get(), csrIA.get(), csrJA.get(),
                     csrValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( rows.releaseDelayed() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::sparseGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "sparseGEVM ( async = " << async << " ) , result = " << alpha << " * x * storage" )

    static LAMAKernel<CSRKernelTrait::sparseGEVM<ValueType> > sparseGEVM;

    ContextPtr loc = this->getContextPtr();
    sparseGEVM.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteAccess<ValueType> wResult( result, loc );

    // result += alpha * x * thisMatrix, can take advantage of row indexes

    IndexType numNonZeroRows = mRowIndexes.size();

    ReadAccess<IndexType> rows( mRowIndexes, loc );

    sparseGEVM[loc]( wResult.get(), alpha, rX.get(),
                     mNumColumns, numNonZeroRows, rows.get(),
                     csrIA.get(), csrJA.get(), csrValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( rows.releaseDelayed() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::normalGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool async ) const
{
    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    // Note: alias &result == &y possible
    //       ReadAccess on y before WriteOnlyAccess on result guarantees valid data

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumValues, csrIA.get(),
                     csrJA.get(), csrValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::normalGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "normalGEVM ( async = " << async << " ) , result = " << alpha << " * x * storage"
                   << " + " << beta << " * y" )

    static LAMAKernel<CSRKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    // Note: alias &result == &y possible
    //       ReadAccess on y before WriteOnlyAccess on result guarantees valid data

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

    normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(),
                     mNumRows, mNumColumns,
                     csrIA.get(), csrJA.get(), csrValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::normalGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

    ValueType beta = 0;

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, NULL, mNumRows, mNumColumns, mNumValues,
                     csrIA.get(), csrJA.get(), csrValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::normalGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "normalGEVM ( async = " << async << " ) , result = " << alpha << " * x * storage" )

    static LAMAKernel<CSRKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );

    ValueType beta = 0;

    normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, NULL,
                     mNumRows, mNumColumns,
                     csrIA.get(), csrJA.get(), csrValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( csrIA.releaseDelayed() );
        syncToken->pushRoutine( csrJA.releaseDelayed() );
        syncToken->pushRoutine( csrValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool  async ) const
{
    SCAI_REGION( "Storage.CSR.gemv" )

    SCAI_LOG_INFO( logger,
                   "GEMV ( async = " << async << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::constants::ZERO || ( mNumValues == 0 ) )
    {
        // so we just have result = beta * y, will be done synchronously

        HArrayUtils::assignScaled( result, beta, y, this->getContextPtr() );

        if ( async )
        {
            return new tasking::NoSyncToken();
        }
        else
        {
            return NULL;
        }
    }

    // check for correct sizes of x

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )

    if ( beta == common::constants::ZERO )
    {
        // take version that does not access y at all (can be undefined or aliased to result)

        return normalGEMV( result, alpha, x, async );
    }

    // y is relevant, so it must have correct size

    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    if ( &result == &y && ( beta == common::constants::ONE ) && ( mRowIndexes.size() > 0 ) )
    {
        // y += A * x,  where only some rows in A are filled, uses more efficient routine

        return sparseGEMV( result, alpha, x, async );
    }
    else
    {
        return normalGEMV( result, alpha, x, beta, y, async );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::gevm(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool  async ) const
{
    SCAI_REGION( "Storage.CSR.gevm" )

    SCAI_LOG_INFO( logger,
                   "GEVM ( async = " << async << " ), result = " << alpha << " * x * A  + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::constants::ZERO || ( mNumValues == 0 ) )
    {
        SCAI_LOG_INFO( logger, "gevm, alpha = 0 : result = " << beta << " * y " )

        // so we just have result = beta * y, will be done synchronously

        HArrayUtils::assignScaled( result, beta, y, this->getContextPtr() );

        if ( async )
        {
            return new tasking::NoSyncToken();
        }
        else
        {
            return NULL;
        }
    }

    // check for correct sizes of x

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )

    if ( beta == common::constants::ZERO )
    {
        // take version that does not access y at all (can be undefined or aliased to result)

        return normalGEVM( result, alpha, x, async );
    }

    // y is relevant, so it must have correct size

    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )

    if ( &result == &y && ( beta == common::constants::ONE ) && ( mRowIndexes.size() > 0 ) )
    {
        // y += x * A,  where only some rows in A are filled, uses more efficient routine

        return sparseGEVM( result, alpha, x, async );
    }
    else
    {
        return normalGEVM( result, alpha, x, beta, y, async );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    bool async = true;

    SyncToken* token = gemv( result, alpha, x, beta, y, async );

    SCAI_ASSERT( token, "NULL token not allowed for asynchronous execution gemv, alpha = " << alpha << ", beta = " << beta )

    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRStorage<ValueType>::vectorTimesMatrixAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    bool async = true;

    SyncToken* token = gevm( result, alpha, x, beta, y, async );

    SCAI_ASSERT( token, "NULL token not allowed for asynchronous execution gevm, alpha = " << alpha << ", beta = " << beta )

    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.CSR.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )

    // matrix must be square

    static LAMAKernel<CSRKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );

    ReadAccess<IndexType> csrIA( mIa, loc );
    ReadAccess<IndexType> csrJA( mJa, loc );
    ReadAccess<ValueType> csrValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );

    // Due to diagonal property there is no advantage by taking row indexes

    SCAI_CONTEXT_ACCESS( loc )

    jacobi[loc]( wSolution.get(), csrIA.get(), csrJA.get(), csrValues.get(),
                 rOldSolution.get(), rRhs.get(), omega, mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.CSR.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumRows() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumColumns() )
    SCAI_ASSERT_DEBUG( localStorage.hasDiagonalProperty(), localStorage << ": has not diagonal property" )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    const CSRStorage<ValueType>* csrLocal;

    if( localStorage.getFormat() == Format::CSR )
    {
        csrLocal = dynamic_cast<const CSRStorage<ValueType>*>( &localStorage );
        SCAI_ASSERT_DEBUG( csrLocal, "could not cast to CSRStorage " << localStorage )
    }
    else
    {
        // either copy localStorage to CSR (not recommended) or
        // just get the diagonal in localValues and set order in localIA
        COMMON_THROWEXCEPTION( "local stroage is not CSR" )
    }

    static LAMAKernel<CSRKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = this->getContextPtr();
    jacobiHalo.getSupportedContext( loc );

    {
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<IndexType> localIA( csrLocal->mIa, loc );
        ReadAccess<ValueType> localValues( csrLocal->mValues, loc );
        ReadAccess<IndexType> haloIA( mIa, loc );
        ReadAccess<IndexType> haloJA( mJa, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

        const IndexType numNonEmptyRows = mRowIndexes.size();

        SCAI_LOG_INFO( logger, "#row indexes = " << numNonEmptyRows )

        if( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );

            SCAI_CONTEXT_ACCESS( loc )

            jacobiHalo[loc]( wSolution.get(), localIA.get(), localValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                             haloRowIndexes.get(), rOldHaloSolution.get(), omega, numNonEmptyRows );
        }
        else
        {
            SCAI_CONTEXT_ACCESS( loc )

            jacobiHalo[loc]( wSolution.get(), localIA.get(), localValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                             NULL, rOldHaloSolution.get(), omega, mNumRows );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.CSR.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    static LAMAKernel<CSRKernelTrait::jacobiHaloWithDiag<ValueType> > jacobiHaloWithDiag;

    ContextPtr loc = this->getContextPtr();
    jacobiHaloWithDiag.getSupportedContext( loc );

    {
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<ValueType> localDiagValues( localDiagonal, loc );
        ReadAccess<IndexType> haloIA( mIa, loc );
        ReadAccess<IndexType> haloJA( mJa, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

        const IndexType numNonEmptyRows = mRowIndexes.size();

        SCAI_LOG_INFO( logger, "#row indexes = " << numNonEmptyRows )

        if( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );

            SCAI_CONTEXT_ACCESS( loc )

            jacobiHaloWithDiag[loc]( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                                     haloRowIndexes.get(), rOldHaloSolution.get(), omega, numNonEmptyRows );
        }
        else
        {
            SCAI_CONTEXT_ACCESS( loc )

            jacobiHaloWithDiag[loc]( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
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
    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", with A = " << a << ", B = " << b )

    SCAI_REGION( "Storage.CSR.plusMatrix" )

    // a and b have to be CSR storages, otherwise create temporaries.

    const CSRStorage<ValueType>* csrA = NULL;
    const CSRStorage<ValueType>* csrB = NULL;

    // Define shared pointers in case we need temporaries

    common::shared_ptr<CSRStorage<ValueType> > tmpA;
    common::shared_ptr<CSRStorage<ValueType> > tmpB;

    if( a.getFormat() == Format::CSR )
    {
        csrA = dynamic_cast<const CSRStorage<ValueType>*>( &a );
        SCAI_ASSERT_DEBUG( csrA, "could not cast to CSRStorage " << a )
    }
    else
    {
        SCAI_UNSUPPORTED( a << ": will be converted to CSR for matrix multiply" )
        tmpA = common::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( a ) );
        csrA = tmpA.get();
    }

    if( b.getFormat() == Format::CSR )
    {
        csrB = dynamic_cast<const CSRStorage<ValueType>*>( &b );
        SCAI_ASSERT_DEBUG( csrB, "could not cast to CSRStorage " << b )
    }
    else
    {
        SCAI_UNSUPPORTED( b << ": will be converted to CSR for matrix multiply" )
        tmpB = common::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( b ) );
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
    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " +* A * B + " << beta << " * C, with " << "A = " << a << ", B = " << b << ", C = " << c )

    SCAI_REGION( "Storage.CSR.timesMatrix" )

    // a and b have to be CSR storages, otherwise create temporaries.

    const CSRStorage<ValueType>* csrA = NULL;
    const CSRStorage<ValueType>* csrB = NULL;
    const CSRStorage<ValueType>* csrC = NULL;

    // Define two shared pointers in case we need temporaries

    common::shared_ptr<CSRStorage<ValueType> > tmpA;
    common::shared_ptr<CSRStorage<ValueType> > tmpB;
    common::shared_ptr<CSRStorage<ValueType> > tmpC;

    if( a.getFormat() == Format::CSR )
    {
        csrA = dynamic_cast<const CSRStorage<ValueType>*>( &a );
        SCAI_ASSERT_DEBUG( csrA, "could not cast to CSRStorage " << a )
    }
    else
    {
        SCAI_UNSUPPORTED( a << ": will be converted to CSR for matrix multiply" )
        tmpA = common::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( a ) );
        csrA = tmpA.get();
    }

    if( b.getFormat() == Format::CSR )
    {
        csrB = dynamic_cast<const CSRStorage<ValueType>*>( &b );
        SCAI_ASSERT_DEBUG( csrB, "could not cast to CSRStorage " << b )
    }
    else
    {
        SCAI_UNSUPPORTED( b << ": will be converted to CSR for matrix multiply" )
        tmpB = common::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( b ) );
        csrB = tmpB.get();
    }

    if( beta != common::constants::ZERO )
    {
        // c temporary needed if not correct format/type or aliased to this

        if( ( c.getFormat() == Format::CSR ) && ( &c != this ) )
        {
            csrC = dynamic_cast<const CSRStorage<ValueType>*>( &c );
            SCAI_ASSERT_DEBUG( csrC, "could not cast to CSRStorage " << c )
        }
        else
        {
            SCAI_UNSUPPORTED( c << ": CSR temporary required for matrix add" )
            tmpC = common::shared_ptr<CSRStorage<ValueType> >( new CSRStorage<ValueType>( c ) );
            csrC = tmpC.get();
        }

    }

    // now we have in any case all arguments as CSR Storage

    ContextPtr loc = Context::getHostPtr();

    if( a.getContextPtr()->getType() == b.getContextPtr()->getType() )
    {
        loc = a.getContextPtr();
    }

    ContextPtr saveContext = getContextPtr();

    CSRStorage<ValueType> tmp1;
    tmp1.matrixTimesMatrixCSR( alpha, *csrA, *csrB, loc );
    tmp1.setContextPtr( loc );

    if( beta != common::constants::ZERO )
    {
        CSRStorage<ValueType> tmp2;
        tmp2.matrixAddMatrixCSR( static_cast<ValueType>(1.0), tmp1, beta, *csrC, loc );
        swap( tmp2 );
    }
    else
    {
        swap( tmp1 );
    }

    this->setContextPtr( saveContext );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::matrixAddMatrixCSR(
    const ValueType alpha,
    const CSRStorage<ValueType>& a,
    const ValueType beta,
    const CSRStorage<ValueType>& b,
    const ContextPtr preferedLoc )
{
    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " * A + " << beta << " * B, with " << "A = " << a << ", B = " << b << ", all are CSR" )

//    // TODO: just temporary, MAKE loc const again!
//    loc = Context::getContextPtr( Context::Host );

    static LAMAKernel<CSRKernelTrait::matrixAddSizes> matrixAddSizes;
    static LAMAKernel<CSRKernelTrait::matrixAdd<ValueType> > matrixAdd;

    ContextPtr loc = preferedLoc;
    matrixAdd.getSupportedContext( loc, matrixAddSizes );

    if( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary

        CSRStorage<ValueType> tmp;

        tmp.matrixAddMatrixCSR( alpha, a, beta, b, loc );

        swap( tmp ); // safe as tmp will be destroyed afterwards

        return;
    }

    SCAI_REGION( "Storage.CSR.addMatrixCSR" )

    allocate( a.getNumRows(), a.getNumColumns() );

    SCAI_ASSERT_EQUAL_ERROR( mNumRows, b.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, b.getNumColumns() )

    mDiagonalProperty = ( mNumRows == mNumColumns );

    {
        ReadAccess<IndexType> aIa( a.getIA(), loc );
        ReadAccess<IndexType> aJa( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );

        ReadAccess<IndexType> bIa( b.getIA(), loc );
        ReadAccess<IndexType> bJa( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );

        // Step 1: compute row sizes of C, build offsets
        SCAI_LOG_DEBUG( logger, "Determing sizes of result matrix C" )

        WriteOnlyAccess<IndexType> cIa( mIa, loc, mNumRows + 1 );

        SCAI_CONTEXT_ACCESS( loc )

        mNumValues = matrixAddSizes[loc] ( cIa.get(), mNumRows, mNumColumns, mDiagonalProperty, aIa.get(), aJa.get(),
                                           bIa.get(), bJa.get() );

        // Step 2: fill in ja, values

        SCAI_LOG_DEBUG( logger, "Compute the sparse values, # = " << mNumValues )

        WriteOnlyAccess<IndexType> cJa( mJa, loc, mNumValues );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValues );

        matrixAdd[loc]( cJa.get(), cValues.get(), cIa.get(), mNumRows, mNumColumns, mDiagonalProperty, alpha, aIa.get(),
                        aJa.get(), aValues.get(), beta, bIa.get(), bJa.get(), bValues.get() );
    }

    SCAI_LOG_DEBUG( logger, *this << ": compress by removing zero elements" )

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
    const ContextPtr preferedLoc )
{
    SCAI_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B, with " << "A = " << a << ", B = " << b << ", all are CSR" << ", Context = " << preferedLoc->getType() )

    // get availabe implementations of needed kernel routines

    static LAMAKernel<CSRKernelTrait::matrixMultiplySizes> matrixMultiplySizes;
    static LAMAKernel<CSRKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    // choose Context where all kernel routines are available

    ContextPtr loc = preferedLoc;
    matrixMultiply.getSupportedContext( loc, matrixMultiplySizes );

    SCAI_ASSERT_ERROR( &a != this, "matrixTimesMatrix: alias of a with this result matrix" )
    SCAI_ASSERT_ERROR( &b != this, "matrixTimesMatrix: alias of b with this result matrix" )

    SCAI_ASSERT_EQUAL_ERROR( a.getNumColumns(), b.getNumRows() )

    IndexType k = a.getNumColumns();

    SCAI_REGION( "Storage.CSR.timesMatrixCSR" )

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

        SCAI_CONTEXT_ACCESS( loc )

        mNumValues = matrixMultiplySizes[loc] ( cIA.get(), mNumRows, mNumColumns, k, mDiagonalProperty, aIA.get(), aJA.get(),
                                                bIA.get(), bJA.get() );

        WriteOnlyAccess<IndexType> cJa( mJa, loc, mNumValues );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValues );

        matrixMultiply[loc]( cIA.get(), cJa.get(), cValues.get(), mNumRows, mNumColumns, k, alpha, mDiagonalProperty,
                             aIA.get(), aJA.get(), aValues.get(), bIA.get(), bJA.get(), bValues.get() );
    }

    // TODO: check this!
//    compress();
    buildRowIndexes();
//    check( "result of matrix x matrix" ); // just verify for a correct matrix
//    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    return utilskernel::HArrayUtils::asum( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    ContextPtr prefLoc = this->getContextPtr();

    ValueType res = HArrayUtils::dotProduct( mValues, mValues, prefLoc );

    return common::Math::sqrt( res );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::maxNorm() const
{
    // no more checks needed here

    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    if( mNumValues == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<UtilKernelTrait::reduce<ValueType> > reduce;

    ContextPtr loc = this->getContextPtr();
    reduce.getSupportedContext( loc );

    ReadAccess<ValueType> csrValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType maxval = reduce[loc]( csrValues.get(), mNumValues, utilskernel::reduction::ABS_MAX );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType CSRStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_REGION( "Storage.CSR.maxDiffNorm" )

    SCAI_ASSERT_EQUAL_ERROR( mNumRows, other.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, other.getNumColumns() )

    SCAI_LOG_INFO( logger, *this << ": maxDiffNorm( " << other << " )" )

    common::shared_ptr<CSRStorage<ValueType> > tmpOtherCSR;

    const CSRStorage<ValueType>* otherCSR;

    if( other.getValueType() == this->getValueType() && ( other.getFormat() == Format::CSR ) )
    {
        otherCSR = dynamic_cast<const CSRStorage<ValueType>*>( &other );
        SCAI_ASSERT_ERROR( otherCSR, other << ": could not cast to " << typeName() )
    }
    else
    {
        SCAI_UNSUPPORTED( other << ": converted to " << typeName() << " for maxDiffNorm" )
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

    SCAI_LOG_INFO( logger, *this << ": maxDiffNormImpl( " << other << " )" )

    if( mNumRows == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    bool sorted = mSortedRows && other.mSortedRows && ( mDiagonalProperty == other.mDiagonalProperty );

    static LAMAKernel<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;

    ContextPtr loc = this->getContextPtr();
    absMaxDiffVal.getSupportedContext( loc );

    ReadAccess<IndexType> csrIA1( mIa, loc );
    ReadAccess<IndexType> csrJA1( mJa, loc );
    ReadAccess<ValueType> csrValues1( mValues, loc );

    ReadAccess<IndexType> csrIA2( other.mIa, loc );
    ReadAccess<IndexType> csrJA2( other.mJa, loc );
    ReadAccess<ValueType> csrValues2( other.mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxDiffVal[loc] ( mNumRows, sorted, csrIA1.get(), csrJA1.get(), csrValues1.get(), csrIA2.get(),
                                            csrJA2.get(), csrValues2.get() );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
CSRStorage<ValueType>* CSRStorage<ValueType>::copy() const
{
    return new CSRStorage<ValueType>( *this );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void CSRStorage<ValueType>::buildSparseRowSizes( HArray<IndexType>& rowSizes ) const
{
    SCAI_LOG_DEBUG( logger, "copy nnz for each row in HArray" );

    WriteOnlyAccess<IndexType> writeRowSizes( rowSizes, mNumRows );
    ReadAccess<IndexType> csrIA( mIa );

    OpenMPCSRUtils::offsets2sizes( writeRowSizes.get(), csrIA.get(), mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CSRStorage<ValueType>::buildSparseRowData(
    HArray<IndexType>& sparseJA,
    HArray<ValueType>& sparseValues ) const
{
    SCAI_LOG_INFO( logger, *this << ": build sparse row data" );

    // for CSR format we can just copy arrays with column indexes and data values

    sparseJA = mJa;
    sparseValues = mValues;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* CSRStorage<ValueType>::create()
{
    return new CSRStorage<ValueType>();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType CSRStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::CSR, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string CSRStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string("CSRStorage<") << common::getScalarType<ValueType>() << std::string(">");
    return s.str();
}

template<typename ValueType>
const char* CSRStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( CSRStorage, SCAI_ARITHMETIC_HOST )

#define CSR_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                  \
    template void CSRStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&, \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );  \
    template  void CSRStorage<ValueType>::setCSRDataSwap(const IndexType, const IndexType, const IndexType,                \
            hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>&,             \
            hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );        \
    template void CSRStorage<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;              \
    template void CSRStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                          \
    template void CSRStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                          \
    template void CSRStorage<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                                \
    template void CSRStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;  \

#define CSR_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, CSR_STORAGE_INST_LVL2, SCAI_ARITHMETIC_HOST )

    SCAI_COMMON_LOOP( CSR_STORAGE_INST_LVL1, SCAI_ARITHMETIC_HOST )

#undef CSR_STORAGE_INST_LVL2
#undef CSR_STORAGE_INST_LVL1


    } /* end namespace lama */

    } /* end namespace scai */
