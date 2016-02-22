/**
 * @file ELLStorage.cpp
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
 * @brief Instantitions for template class ELLStorage.
 * @author Lauretta Schubert
 * @date 25.05.2011
 */

// hpp
#include <scai/lama/storage/ELLStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/ELLKernelTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/hmemo.hpp>

#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/exception/UnsupportedException.hpp>
#include <scai/common/preprocessor.hpp>

namespace scai
{

using common::shared_ptr;
using common::unique_ptr;
using tasking::SyncToken;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;
using utilskernel::LArray;

using sparsekernel::ELLKernelTrait;
using sparsekernel::CSRKernelTrait;

using namespace tasking;
using namespace hmemo;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ELLStorage<ValueType>::logger, "MatrixStorage.ELLStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const common::context::ContextType con /* = Context::Host */)

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>( numRows, numColumns ), mNumValuesPerRow( 0 )
{
    // TODO in other formats the last parameter is "const ContextPtr loc"

    ContextPtr loc = Context::getContextPtr( con );

    setContextPtr( loc );

    // Initialization requires correct values for the IA array with 0

    static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;

    WriteOnlyAccess<IndexType> ellSizes( mIA, loc, mNumRows );

    SCAI_CONTEXT_ACCESS( loc )

    setVal[loc]( ellSizes.get(), mNumRows, 0, common::reduction::COPY );

    SCAI_LOG_DEBUG( logger, "ELLStorage for matrix " << mNumRows << " x " << mNumColumns << ", no elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage()

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>( 0, 0 ), mNumValuesPerRow( 0 )
{
    SCAI_LOG_DEBUG( logger, "ELLStorage, default constructor for zero matrix." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRows,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values )

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>()
{
    SCAI_LOG_INFO( logger, "constructor with ELL data array" )

    setELLData( numRows, numColumns, numValuesPerRows, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( const ELLStorage<ValueType>& other )

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>( 0, 0 )
{
    SCAI_LOG_INFO( logger, "constructor # other = " << other )

    // call the assignment operator

    operator=( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>& ELLStorage<ValueType>::operator=( const ELLStorage<ValueType>& other )
{
    // nothing to do for self assignments

    if( &other == this )
    {
        return *this;
    }

    // For optimization: fillValues, resetDiagonalProperty, check is not really needed

    setELLData( other.mNumRows, other.mNumColumns, other.mNumValuesPerRow, other.mIA, other.mJA, other.mValues );

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>& ELLStorage<ValueType>::operator=( const _MatrixStorage& other )
{
    assign( other ); // calls virtual method of MatrixStorage

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::print() const
{
    SCAI_LOG_INFO( logger, "print" )

    using std::cout;
    using std::endl;

    cout << "ELLStorage " << mNumRows << " x " << mNumColumns << ", #values = " << getNumValues() << endl;

    ReadAccess<IndexType> ia( mIA );
    ReadAccess<IndexType> ja( mJA );
    ReadAccess<ValueType> values( mValues );

    for( IndexType i = 0; i < mNumRows; i++ )
    {
        cout << "Row " << i << " ( " << ia[i] << " entries ) :";

        for( IndexType jj = 0; jj < ia[i]; ++jj )
        {
            IndexType pos = jj * mNumRows + i;
            cout << " " << ja[pos] << ":" << values[pos];
        }

        cout << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorageFormat ELLStorage<ValueType>::getFormat() const
{
    return Format::ELL;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType ELLStorage<ValueType>::getNumValues() const
{
    SCAI_LOG_INFO( logger, "getNumValues" )

    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;

    ContextPtr loc = reduce.getValidContext( this->getContextPtr() );

    ReadAccess<IndexType> ia( mIA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType numValues = reduce[ loc ]( ia.get(), mNumRows, common::reduction::ADD );

    return numValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::purge()
{
    SCAI_LOG_INFO( logger, "purge" )

    mNumColumns = 0;
    mNumRows = 0;
    mNumValuesPerRow = 0;

    mIA.purge();
    mJA.purge();
    mValues.purge();
    mRowIndexes.purge();

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_INFO( logger, "set identity # size = " << size )

    mNumRows = size;
    mNumColumns = size;
    mNumValuesPerRow = 1;

    {
        static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;

        ContextPtr loc = setVal.getValidContext( this->getContextPtr() );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> ia( mIA, loc, mNumRows );

        setVal[loc]( ia.get(), mNumRows, 1, common::reduction::COPY );
    }

    {
        static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;

        ContextPtr loc = setOrder.getValidContext( this->getContextPtr() );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> ja( mJA, loc, mNumRows );

        setOrder[loc]( ja.get(), mNumRows );
    }

    // extra block caused by differnt types of setVal()
    {
        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

        ContextPtr loc = setVal.getValidContext( this->getContextPtr() );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<ValueType> data( mValues, loc, mNumRows );

        setVal[loc]( data.get(), mNumRows, ValueType( 1 ), common::reduction::COPY );
    }

    mDiagonalProperty = true;

    SCAI_LOG_INFO( logger, *this << " is identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool ELLStorage<ValueType>::checkDiagonalProperty() const
{
    SCAI_LOG_INFO( logger, "checkDiagonalProperty" )

    IndexType numDiagonals = std::min( mNumRows, mNumColumns );

    bool diagonalProperty = true;

    if( numDiagonals == 0 )
    {
        // diagonal property is given for zero-sized matrices

        diagonalProperty = true;
    }
    else if( mNumValuesPerRow < 1 )
    {
        // no elements, so certainly it does not have diagonl property

        diagonalProperty = false;
    }
    else
    {
        static LAMAKernel<ELLKernelTrait::hasDiagonalProperty> ellHasDiagonalProperty;

        // check it where the JA array has a valid copy

        ContextPtr loc = ellHasDiagonalProperty.getValidContext( mJA.getValidContext() );

        ReadAccess<IndexType> ja( mJA, loc );

        SCAI_CONTEXT_ACCESS( loc )

        diagonalProperty = ellHasDiagonalProperty[loc]( numDiagonals, ja.get() );
    }

    SCAI_LOG_INFO( logger, *this << ": checkDiagonalProperty = " << diagonalProperty )

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::clear()
{
    SCAI_LOG_INFO( logger, "clear" )

    mNumRows = 0;
    mNumColumns = 0;
    mNumValuesPerRow = 0;

    mIA.clear();
    mJA.clear();
    mValues.clear();
    mRowIndexes.clear();

    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<OtherValueType>* values,
    const ContextPtr context ) const
{
    SCAI_REGION( "Storage.ELL->CSR" )

    SCAI_LOG_INFO( logger,
                   "buildCSR<" << common::getScalarType<OtherValueType>() << ">" 
                    << " from ELL<" << common::getScalarType<ValueType>() << ">" 
                    << " on " << *context << " ( preferred )" )

    // step 1 : compute IA offsets

    IndexType numValues = 0;

    {
        static LAMAKernel<UtilKernelTrait::set<IndexType, IndexType> > set;
        static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;

        const ContextPtr loc = set.getValidContext( sizes2offsets, set, context );
  
        ReadAccess<IndexType> ellSizes( mIA, loc );
        WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );

        SCAI_CONTEXT_ACCESS( loc )

        // just copy the size array mIA

        set[loc]( csrIA.get(), ellSizes.get(), mNumRows, common::reduction::COPY );

        if( ja == NULL || values == NULL )
        {
            csrIA.resize( mNumRows );
            return;
        }

        numValues = sizes2offsets[loc]( csrIA.get(), mNumRows );
    }

    static LAMAKernel<ELLKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;

    const ContextPtr loc = getCSRValues.getValidContext( context );

    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<IndexType> csrIA( ia, loc );
    ReadAccess<IndexType> ellSizes( mIA, loc );

    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<OtherValueType> csrValues( *values, loc, numValues );

    SCAI_CONTEXT_ACCESS( loc )

    getCSRValues[loc]( csrJA.get(), csrValues.get(), csrIA.get(), mNumRows, mNumValuesPerRow, 
                       ellSizes.get(), ellJA.get(), ellValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr context )
{
    SCAI_REGION( "Storage.ELL<-CSR" )

    SCAI_LOG_INFO( logger,
                   "set CSR data on " << *context << ": numRows = " << numRows << ", numColumns = " << numColumns 
                   << ", numValues = " << numValues << ", compress threshold = " << mCompressThreshold )

    if( numRows == 0 )
    {
        // just allocate will clear member arrays

        allocate( numRows, numColumns );

        return;
    }

    _MatrixStorage::setDimension( numRows, numColumns );

    // Get function pointers for needed routines at the LAMA interface

    static LAMAKernel<CSRKernelTrait::offsets2sizes > offsets2sizes;
    static LAMAKernel<ELLKernelTrait::hasDiagonalProperty > hasDiagonalProperty;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<ELLKernelTrait::setCSRValues<ValueType, OtherValueType> > setCSRValues;

    ContextPtr loc = offsets2sizes.getValidContext( hasDiagonalProperty, setCSRValues, context );

    // build array with non-zero values per row

    {
        ReadAccess<IndexType> csrIA( ia, loc );
        WriteOnlyAccess<IndexType> ellSizes( mIA, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )
        offsets2sizes[ loc ]( ellSizes.get(), csrIA.get(), mNumRows );
    }

    // determine the maximal number of non-zero in one row

    {
        ReadAccess<IndexType> ellSizes( mIA, loc );
        SCAI_CONTEXT_ACCESS( loc )
        mNumValuesPerRow = reduce[loc]( ellSizes.get(), mNumRows, common::reduction::MAX );
    }

    SCAI_LOG_INFO( logger, "setCSRData, #values/row = " << mNumValuesPerRow )

    //  Now we know the size of the ja and values arrays for the ELL format

    const IndexType dataSize = mNumValuesPerRow * mNumRows;

    if( mNumRows > 200 && mNumValuesPerRow > 0 )
    {
        // make this check only on larger matrices, dataSize must not be equal 0

        double fillRate = double( numValues ) / double( dataSize );

        if( fillRate < 0.5 )
        {
            SCAI_LOG_WARN( logger,
                           *this << ": fill rate = " << fillRate << " ( " << numValues << " non-zero values ), consider using JDS" )
        }
    }

    {
        // now fill the matrix values and column indexes

        ReadAccess<IndexType> csrIA( ia, loc );
        ReadAccess<IndexType> csrJA( ja, loc );
        ReadAccess<OtherValueType> csrValues( values, loc );

        ReadAccess<IndexType> ellIA( mIA, loc );

        WriteOnlyAccess<IndexType> ellJA( mJA, loc, dataSize );
        WriteOnlyAccess<ValueType> ellValues( mValues, loc, dataSize );

        SCAI_LOG_DEBUG( logger, "convert CSR -> ELL, ellSize = " << dataSize )

        SCAI_CONTEXT_ACCESS( loc )

        setCSRValues[loc]( ellJA.get(), ellValues.get(), ellIA.get(), 
                           mNumRows, mNumValuesPerRow, 
                           csrIA.get(), csrJA.get(), csrValues.get() );

        SCAI_LOG_DEBUG( logger, " size = " <<ellJA.size() )

        IndexType numDiagonals = std::min( mNumRows, mNumColumns );

        if( numDiagonals == 0 )
        {
            mDiagonalProperty = true;
        }
        else if( numValues == 0 )
        {
            mDiagonalProperty = false;
        }
        else
        {
            SCAI_CONTEXT_ACCESS( loc )
            mDiagonalProperty = hasDiagonalProperty[loc]( numDiagonals, ellJA.get() );
        }
    }

    if( numRows == numColumns && !mDiagonalProperty )
    {
        SCAI_LOG_INFO( logger, *this << ": square matrix has not diagonal property" )
    }

    buildRowIndexes( loc );

    SCAI_LOG_DEBUG( logger, "convert CSR -> ELL done: " << *this )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void ELLStorage<ValueType>::setELLData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    SCAI_ASSERT_EQUAL_ERROR( numRows, ia.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows * numValuesPerRow, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows * numValuesPerRow, values.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValuesPerRow = numValuesPerRow;

    ContextPtr loc = getContextPtr();

    HArrayUtils::setImpl( mIA, ia, common::reduction::COPY, loc );
    HArrayUtils::setImpl( mJA, ja, common::reduction::COPY, loc );

    HArrayUtils::set( mValues, values, common::reduction::COPY, loc );  // also type conversion

    // fill up my arrays ja and values to make matrix-multiplication fast

    {
        static LAMAKernel<ELLKernelTrait::fillELLValues<ValueType> > fillELLValues;

        ContextPtr loc = fillELLValues.getValidContext( this->getContextPtr() );

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> ellIA( mIA, loc );
        WriteAccess<IndexType> ellJA( mJA, loc );
        WriteAccess<ValueType> ellValues( mValues, loc );

        fillELLValues[loc]( ellJA.get(), ellValues.get(), ellIA.get(), mNumRows, mNumValuesPerRow );
    }

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode

#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "ELLStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif

    this->resetDiagonalProperty();

    SCAI_LOG_INFO( logger, *this << ": set ELLPACK by arrays ia, ja, values" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType ELLStorage<ValueType>::getNumValuesPerRow() const
{
    SCAI_LOG_INFO( logger, "getNumValuesPerRow" )
    return mNumValuesPerRow;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    SCAI_LOG_INFO( logger, "setDiagonalImpl # value = " << value )

    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;

    ContextPtr loc = setVal.getValidContext( this->getContextPtr() );

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    SCAI_CONTEXT_ACCESS( loc )

    WriteAccess<ValueType> wValues( mValues, loc );

    setVal[ loc ]( wValues.get(), numDiagonalElements, value, common::reduction::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::setDiagonalImpl( const HArray<OtherType>& diagonal )
{
    SCAI_LOG_INFO( logger, "setDiagonalImpl # diagonal = " << diagonal )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    static LAMAKernel<UtilKernelTrait::set<ValueType, OtherType> > set;

    ContextPtr loc = set.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    // ELL format with diagonal property: diagonal is just the first column in mValues

    set[ loc ]( wValues.get(), rDiagonal.get(), numDiagonalElements, common::reduction::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_LOG_TRACE( logger, "getRowImpl # row = " << row << ", i = " << i )

    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    static LAMAKernel<ELLKernelTrait::getRow<ValueType, OtherType> > getRow;

    ContextPtr loc = getRow.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<OtherType> wRow( row, loc, mNumColumns );

    const ReadAccess<IndexType> rIa( mIA, loc );
    const ReadAccess<IndexType> rJa( mJA, loc );
    const ReadAccess<ValueType> rValues( mValues, loc );

    getRow[loc]( wRow.get(), i, mNumRows, mNumColumns, mNumValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::getDiagonalImpl( HArray<OtherType>& diagonal ) const
{
    SCAI_LOG_INFO( logger, "getDiagonalImpl # diagonal = " << diagonal )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    // OtherType is output type, so use it as first template argument

    static LAMAKernel<UtilKernelTrait::set<OtherType, ValueType> > set;

    ContextPtr loc = set.getValidContext( this->getContextPtr() );

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );

    // ELL format with diagonal property: diagonal is just the first column in mValues

    SCAI_CONTEXT_ACCESS( loc )

    set[loc]( wDiagonal.get(), rValues.get(), numDiagonalElements, common::reduction::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::scaleImpl( const ValueType value )
{
    SCAI_LOG_INFO( logger, "scaleImpl # value = " << value )

    HArrayUtils::scale( mValues, value, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::conj()
{
    HArrayUtils::conj( mValues, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::scaleImpl( const HArray<OtherValueType>& values )
{
    SCAI_LOG_INFO( logger, "scaleImpl # values = " << values )

    static LAMAKernel<ELLKernelTrait::scaleValue<ValueType, OtherValueType> > scaleValue;

    ContextPtr loc = scaleValue.getValidContext( this->getContextPtr() );

    ReadAccess<OtherValueType> rValues( values, loc );
    ReadAccess<IndexType> rIa( mIA, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    scaleValue[loc]( mNumRows, mNumValuesPerRow, rIa.get(), wValues.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LArray<IndexType>& ELLStorage<ValueType>::getIA() const
{
    return mIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LArray<IndexType>& ELLStorage<ValueType>::getJA() const
{
    return mJA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LArray<ValueType>& ELLStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::~ELLStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~ELLStorage for matrix " << mNumRows << " x " << mNumColumns << ", # nnr = " << mNumValuesPerRow )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::check( const char* msg ) const
{
    SCAI_LOG_INFO( logger, "check # msg = " << msg )

    SCAI_ASSERT_EQUAL_ERROR( mNumRows, mIA.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValuesPerRow * mNumRows, mJA.size() )
    SCAI_ASSERT_EQUAL_ERROR( mJA.size(), mValues.size() )

    static LAMAKernel<ELLKernelTrait::check> check;

    ContextPtr loc = check.getValidContext( this->getContextPtr() );

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    check[loc]( mNumRows, mNumValuesPerRow, mNumColumns, rIa.get(), rJa.get(), msg );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate ELL sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    mNumRows = numRows;
    mNumColumns = numColumns;

    SCAI_LOG_DEBUG( logger, "resize mIA, mNumRows = " << mNumRows )

    {
        // Intialize array mIA with 0

        static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;

        ContextPtr loc = setVal.getValidContext( getContextPtr() );

        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> ia( mIA, loc, mNumRows );

        setVal[ loc ]( ia.get(), mNumRows, 0, common::reduction::COPY );
    }

    mDiagonalProperty = checkDiagonalProperty();

    SCAI_LOG_DEBUG( logger, "ready allocate" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ELLStorage<" << common::getScalarType<ValueType>() 
           << ">( size = " << mNumRows << " x " << mNumColumns
           << ", nnr = " << mNumValuesPerRow << ", threshold = " << mCompressThreshold << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )
    SCAI_LOG_TRACE( logger, "sizes: ia = " << mIA.size() << ", ja = " << mJA.size() << ", data = " << mValues.size() )

    static LAMAKernel<ELLKernelTrait::getValue<ValueType> > getValue;

    ContextPtr loc = getValue.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    const ReadAccess<IndexType> rIa( mIA, loc );
    const ReadAccess<IndexType> rJa( mJA, loc );
    const ReadAccess<ValueType> rValues( mValues, loc );

    return getValue[loc]( i, j, mNumRows, mNumValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::prefetch( const ContextPtr context ) const
{
    SCAI_LOG_INFO( logger, "prefetch # context " << context )
    SCAI_LOG_DEBUG( logger, "Starting prefetch of "<<*this<<" to "<<context )

    mRowIndexes.prefetch( context );
    mIA.prefetch( context );
    mJA.prefetch( context );
    mValues.prefetch( context );

    SCAI_LOG_DEBUG( logger, "Finished prefetch of " << *this << " to " << context )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::wait() const
{
    SCAI_LOG_INFO( logger, "wait" )

    mRowIndexes.wait();
    mIA.wait();
    mJA.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::buildRowIndexes( const ContextPtr context )
{
    SCAI_LOG_INFO( logger, "buildRowIndexes # loc = " << context )

    mRowIndexes.clear();

    if ( mNumRows == 0 )
    {
        return;
    }

    // Get function pointers for needed kernel routines

    static LAMAKernel<ELLKernelTrait::countNonEmptyRowsBySizes> countNonEmptyRowsBySizes;
    static LAMAKernel<ELLKernelTrait::setNonEmptyRowsBySizes> setNonEmptyRowsBySizes;

    // choose location where both routines are available

    ContextPtr loc = countNonEmptyRowsBySizes.getValidContext( setNonEmptyRowsBySizes, context );

    ReadAccess<IndexType> ellIA( mIA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType nonZeroRows = countNonEmptyRowsBySizes[loc]( ellIA.get(), mNumRows );

    float usage = float( nonZeroRows ) / float( mNumRows );

    if( usage >= mCompressThreshold )
    {
        SCAI_LOG_INFO( logger,
                       "ELLStorage: do not build row indexes, usage = " << usage << " >= " << mCompressThreshold << " ( threshold )" )
        return;
    }

    WriteOnlyAccess<IndexType> rowIndexes( mRowIndexes, loc, nonZeroRows );

    setNonEmptyRowsBySizes[loc]( rowIndexes.get(), nonZeroRows, ellIA.get(), mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::compress( const ValueType eps /* = 0.0 */)
{
    SCAI_LOG_INFO( logger, "compress: eps = " << eps )

    static LAMAKernel<ELLKernelTrait::compressIA<ValueType> > compressIA;
    static LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;

    ContextPtr loc = compressIA.getValidContext( compressValues, this->getContextPtr() );

    ReadAccess<IndexType> IA( mIA, loc );
    ReadAccess<IndexType> JA( mJA, loc );
    ReadAccess<ValueType> values( mValues, loc );

    // 1. Step: Check for 0 elements and write new IA array
    LArray<IndexType> newIAArray;
    WriteOnlyAccess<IndexType> newIA( newIAArray, loc, mNumRows );

    compressIA[loc]( IA.get(), JA.get(), values.get(), mNumRows, mNumValuesPerRow, eps, newIA.get() );

    // 2. Step: compute length of longest row
    IndexType newNumValuesPerRow = reduce[ loc ]( IA.get(), mNumRows, common::reduction::MAX );

    // Do further steps, if new array could be smaller
    if( newNumValuesPerRow < mNumValuesPerRow )
    {
        // 3. Step: Allocate new JA and Values array
        LArray<ValueType> newValuesArray;
        LArray<IndexType> newJAArray;
        WriteOnlyAccess<ValueType> newValues( newValuesArray, loc, mNumRows * newNumValuesPerRow );
        WriteOnlyAccess<IndexType> newJA( newJAArray, loc, mNumRows * newNumValuesPerRow );

        // 4. Step: Compute new JA and Values array
        compressValues[loc]( IA.get(), JA.get(), values.get(), mNumRows, mNumValuesPerRow, eps, newNumValuesPerRow,
                             newJA.get(), newValues.get() );

        mJA.swap( newJAArray );
        mValues.swap( newValuesArray );
        mNumValuesPerRow = newNumValuesPerRow;
    }
}

template<typename ValueType>
void ELLStorage<ValueType>::swap( ELLStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "swap # other = " << other )

    std::swap( mNumValuesPerRow, other.mNumValuesPerRow );
    mIA.swap( other.mIA );
    mJA.swap( other.mJA );
    mValues.swap( other.mValues );

    MatrixStorage<ValueType>::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t ELLStorage<ValueType>::getMemoryUsageImpl() const
{
    SCAI_LOG_INFO( logger, "getMemoryUsageImpl" )

    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType) * mIA.size();
    memoryUsage += sizeof(IndexType) * mJA.size();
    memoryUsage += sizeof(ValueType) * mValues.size();

    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::matrixTimesVector(
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
void ELLStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.ELL.VectorTimesMatrix" )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    static LAMAKernel<ELLKernelTrait::sparseGEVM<ValueType> > sparseGEVM;
    static LAMAKernel<ELLKernelTrait::normalGEVM<ValueType> > normalGEVM;

    ContextPtr loc = sparseGEVM.getValidContext( normalGEVM, this->getContextPtr() );

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )

    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( &result == &y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        if( mRowIndexes.size() > 0 && ( beta == scai::common::constants::ONE ) )
        {
            // y += alpha * thisMatrix * x, can take advantage of row indexes

            IndexType numNonZeroRows = mRowIndexes.size();
            ReadAccess<IndexType> rows( mRowIndexes, loc );

            SCAI_CONTEXT_ACCESS( loc )
            sparseGEVM[loc]( wResult.get(), alpha, rX.get(), mNumRows, mNumColumns, mNumValuesPerRow, numNonZeroRows,
                             rows.get(), ellSizes.get(), ellJA.get(), ellValues.get() );
        }
        else
        {
            // we assume that normalGEVM can deal with the alias of result, y

            SCAI_CONTEXT_ACCESS( loc )
            normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumValuesPerRow,
                             ellSizes.get(), ellJA.get(), ellValues.get() );
        }
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )
        normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumValuesPerRow,
                         ellSizes.get(), ellJA.get(), ellValues.get() );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::matrixTimesVectorAsync(
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
SyncToken* ELLStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool  async ) const
{
    SCAI_REGION( "Storage.ELL.gemv" )

    SCAI_LOG_INFO( logger,
                   "GEMV ( async = " << async << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::constants::ZERO || ( mNumValuesPerRow == 0 ) )
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
SyncToken* ELLStorage<ValueType>::normalGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = normalGEMV.getValidContext( normalGEMV, this->getContextPtr() );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
        syncToken->setCurrent();
    }

    SCAI_CONTEXT_ACCESS( loc )

    // Note: alias &result == &y possible
    //       ReadAccess on y before WriteOnlyAccess on result guarantees valid data

    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumValuesPerRow,
                     ellIA.get(), ellJA.get(), ellValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( ellIA.releaseDelayed() );
        syncToken->pushRoutine( ellJA.releaseDelayed() );
        syncToken->pushRoutine( ellValues.releaseDelayed() );
        syncToken->unsetCurrent();
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::normalGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = normalGEMV.getValidContext( normalGEMV, this->getContextPtr() );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
        syncToken->setCurrent();
    }

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), 0, NULL, mNumRows, mNumValuesPerRow,
                     ellIA.get(), ellJA.get(), ellValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( ellIA.releaseDelayed() );
        syncToken->pushRoutine( ellJA.releaseDelayed() );
        syncToken->pushRoutine( ellValues.releaseDelayed() );
        syncToken->unsetCurrent();
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::sparseGEMV(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = sparseGEMV.getValidContext( sparseGEMV, this->getContextPtr() );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
        syncToken->setCurrent();
    }

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteAccess<ValueType> wResult( result, loc );

    // result += alpha * thisMatrix * x, can take advantage of row indexes

    IndexType numNonZeroRows = mRowIndexes.size();

    ReadAccess<IndexType> rRowIndexes( mRowIndexes, loc );

    sparseGEMV[loc]( wResult.get(), alpha, rX.get(), mNumRows, mNumValuesPerRow, numNonZeroRows,
                     rRowIndexes.get(), ellIA.get(), ellJA.get(), ellValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( rRowIndexes.releaseDelayed() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( ellIA.releaseDelayed() );
        syncToken->pushRoutine( ellJA.releaseDelayed() );
        syncToken->pushRoutine( ellValues.releaseDelayed() );
        syncToken->unsetCurrent();
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::vectorTimesMatrixAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    SCAI_REGION( "Storage.ELL.vectorTimesMatrixAsync" )

    static LAMAKernel<ELLKernelTrait::sparseGEVM<ValueType> > sparseGEVM;
    static LAMAKernel<ELLKernelTrait::normalGEVM<ValueType> > normalGEVM;

    // default location is context of this storage

    ContextPtr loc = normalGEVM.getValidContext( sparseGEVM, this->getContextPtr() );

    // Note: checks will be done by asynchronous task in any case
    //       and exception in tasks are handled correctly

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on " << *loc )

    if( loc->getType() == common::context::Host )
    {
        // execution as separate thread

        void (ELLStorage::*pf)(
            HArray<ValueType>&,
            const ValueType,
            const HArray<ValueType>&,
            const ValueType,
            const HArray<ValueType>& ) const

            = &ELLStorage<ValueType>::vectorTimesMatrix;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrixAsync on Host by own thread" )

        return new tasking::TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( result.size(), mNumColumns )

    if ( ( beta != scai::common::constants::ZERO ) && ( &result != &y ) )
    {
        SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    }

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    SCAI_ASYNCHRONOUS( *syncToken )

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( &result == &y )
    {
        // only write access for y, no read access for result

        WriteAccess<ValueType> wResult( result, loc );

        if( mRowIndexes.size() > 0 && ( beta == scai::common::constants::ONE ) )
        {
            // y += alpha * thisMatrix * x, can take advantage of row indexes

            IndexType numNonZeroRows = mRowIndexes.size();

            ReadAccess<IndexType> rows( mRowIndexes, loc );

            SCAI_CONTEXT_ACCESS( loc )

            sparseGEVM[loc]( wResult.get(), alpha, rX.get(), mNumRows, mNumColumns, mNumValuesPerRow, numNonZeroRows,
                             rows.get(), ellSizes.get(), ellJA.get(), ellValues.get() );

            syncToken->pushRoutine( rows.releaseDelayed() );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            SCAI_CONTEXT_ACCESS( loc )

            normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumColumns, mNumValuesPerRow,
                             ellSizes.get(), ellJA.get(), ellValues.get() );
        }

        syncToken->pushRoutine( wResult.releaseDelayed() );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumValuesPerRow,
                         ellSizes.get(), ellJA.get(), ellValues.get() );

        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
    }

    syncToken->pushRoutine( ellSizes.releaseDelayed() );
    syncToken->pushRoutine( ellJA.releaseDelayed() );
    syncToken->pushRoutine( ellValues.releaseDelayed() );
    syncToken->pushRoutine( rX.releaseDelayed() );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.ELL.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )

    // matrix must be square

    static LAMAKernel<ELLKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = jacobi.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    // make all needed data available at loc

    WriteAccess<ValueType> wSolution( solution, loc );
    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    jacobi[loc] ( wSolution.get(), mNumRows, mNumValuesPerRow, ellSizes.get(), ellJA.get(), ellValues.get(),
                  rOldSolution.get(), rRhs.get(), omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::jacobiIterateAsync(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.ELL.jacobiIterateAsync" )

    static LAMAKernel<ELLKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = jacobi.getValidContext( this->getContextPtr() );

    if ( loc->getType() == common::context::Host )
    {
        // used later in OpenMP to generate a TaskSyncToken

        void (ELLStorage::*jb)(
            HArray<ValueType>&,
            const HArray<ValueType>&,
            const HArray<ValueType>&,
            const ValueType omega ) const

            = &ELLStorage<ValueType>::jacobiIterate;

        using scai::common::bind;
        using scai::common::cref;
        using scai::common::ref;

        return new tasking::TaskSyncToken( bind( jb, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    // matrix must be square, solution vectors must have right size

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    SCAI_ASYNCHRONOUS( *syncToken )

    // make all needed data available at loc

    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    WriteAccess<ValueType> wSolution(  solution, loc );

    SCAI_CONTEXT_ACCESS( loc )

    jacobi[loc]( wSolution.get(), mNumRows, mNumValuesPerRow, ellSizes.get(), ellJA.get(), ellValues.get(),
                 rOldSolution.get(), rRhs.get(), omega );

    syncToken->pushRoutine( rRhs.releaseDelayed() );
    syncToken->pushRoutine( rOldSolution.releaseDelayed() );
    syncToken->pushRoutine( ellValues.releaseDelayed() );
    syncToken->pushRoutine( ellJA.releaseDelayed() );
    syncToken->pushRoutine( ellSizes.releaseDelayed() );
    syncToken->pushRoutine( wSolution.releaseDelayed() );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const HArray<ValueType>& haloOldSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.ELL.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, "HOST: Jacobi iteration on halo matrix data." )

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumRows() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumColumns() )
    SCAI_ASSERT_DEBUG( localStorage.hasDiagonalProperty(), localStorage << ": has not diagonal property" )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, haloOldSolution.size() )

    const HArray<ValueType>* localDiagonal;

    // might be we need a temporary LAMA array for the local diagonal

    common::shared_ptr<HArray<ValueType> > tmpLocalDiagonal;

    if( localStorage.getFormat() == Format::ELL )
    {
        const ELLStorage<ValueType>* ellLocal;

        ellLocal = dynamic_cast<const ELLStorage<ValueType>*>( &localStorage );
        SCAI_ASSERT_DEBUG( ellLocal, "could not cast to ELLStorage " << localStorage )
        localDiagonal = &( ellLocal->mValues );
    }
    else
    {
        // make a temporary for the diagonal and get it from local storage

        SCAI_LOG_WARN( logger, "local stroage is not ELL, temorary needed for diagonal" )

        tmpLocalDiagonal = common::shared_ptr<HArray<ValueType> >( new HArray<ValueType>() );
        localStorage.getDiagonal( *tmpLocalDiagonal );
        localDiagonal = tmpLocalDiagonal.get();

        // Note: tmpLocalDiagonal will be freed at end of routine
    }

    jacobiIterateHalo( localSolution, *localDiagonal, haloOldSolution, omega );

}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& haloOldSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.ELL.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, "HOST: Jacobi iteration on halo matrix data." )

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, haloOldSolution.size() )

    static LAMAKernel<ELLKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = jacobiHalo.getValidContext( this->getContextPtr() );

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<ValueType> rLocalDiagonal( localDiagonal, loc );
        ReadAccess<IndexType> haloIA( mIA, loc );
        ReadAccess<IndexType> haloJA( mJA, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( haloOldSolution, loc );

        const IndexType numNonEmptyRows = mRowIndexes.size();

        if( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );

            jacobiHalo[loc]( wSolution.get(), mNumRows, rLocalDiagonal.get(), mNumValuesPerRow, haloIA.get(), haloJA.get(),
                             haloValues.get(), haloRowIndexes.get(), numNonEmptyRows, rOldHaloSolution.get(), omega );
        }
        else
        {
            // no row indexes available, computation is done over all rows

            const IndexType numNonEmptyRows = mNumRows;

            jacobiHalo[loc]( wSolution.get(), mNumRows, rLocalDiagonal.get(), mNumValuesPerRow, haloIA.get(), haloJA.get(),
                             haloValues.get(), NULL, numNonEmptyRows, rOldHaloSolution.get(), omega );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::l1Norm() const
{
	SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    if( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;

    ContextPtr loc = asum.getValidContext( this->getContextPtr() );

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return asum[loc]( mValues.size(), data.get(), 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::l2Norm() const
{
	SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    if( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;

    ContextPtr loc = dot.getValidContext( this->getContextPtr() );

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc );

	return common::Math::sqrt(dot[loc]( mValues.size(), data.get(), 1, data.get(), 1 ));
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    if( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return static_cast<ValueType>(0.0);
    }

    static LAMAKernel<ELLKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = absMaxVal.getValidContext( this->getContextPtr() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );

    ValueType maxval = absMaxVal[loc]( mNumRows, mNumValuesPerRow, ellIA.get(), ellValues.get() );

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& c )
{
    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " * A * B + " << beta << " * C, with " << "A = " << a << ", B = " << b << ", C = " << c )

    const ELLStorage<ValueType>* ellA = NULL;
    const ELLStorage<ValueType>* ellB = NULL;
    const ELLStorage<ValueType>* ellC = NULL;

    //    common::shared_ptr<CSRStorage<ValueType> > tmpA;
    //    common::shared_ptr<CSRStorage<ValueType> > tmpB;
    common::shared_ptr<ELLStorage<ValueType> > tmpC;

    if( a.getFormat() == Format::ELL )
    {
        ellA = dynamic_cast<const ELLStorage<ValueType>*>( &a );
        SCAI_ASSERT_DEBUG( ellA, "could not cast to ELLStorage " << a )
    }
    else
    {
        SCAI_LOG_ERROR( logger, a << ": a not ELL format" )
    }

    if( b.getFormat() == Format::ELL )
    {
        ellB = dynamic_cast<const ELLStorage<ValueType>*>( &b );
        SCAI_ASSERT_DEBUG( ellB, "could not cast to ELLStorage " << b )
    }
    else
    {
        SCAI_UNSUPPORTED( b << ": b not ELL format" )
    }

    if( ellA == NULL || ellB == NULL )
    {
        // input matrices not ELL format, so try via CSR

        MatrixStorage<ValueType>::matrixTimesMatrix( alpha, a, b, beta, c );
        return;
    }

    if( beta != scai::common::constants::ZERO )
    {
        if( ( c.getFormat() == Format::ELL ) && ( &c != this ) )
        {
            ellC = dynamic_cast<const ELLStorage<ValueType>*>( &c );
            SCAI_ASSERT_DEBUG( ellC, "could not cast to ELLStorage " << c )
        }
        else
        {
            SCAI_UNSUPPORTED( c << ": ELL temporary required for matrix add" )
            tmpC = common::shared_ptr<ELLStorage<ValueType> >( new ELLStorage<ValueType>( c ) );
            ellC = tmpC.get();
        }

    }

    ELLStorage<ValueType> tmp;
    tmp.matrixTimesMatrixELL( alpha, *ellA, *ellB );

    if( beta != scai::common::constants::ZERO )
    {
        ELLStorage<ValueType> tmp1;
        tmp1.matrixAddMatrixELL( static_cast<ValueType>(1.0), tmp, beta, *ellC );
        swap( tmp1 );
    }
    else
    {
        swap( tmp );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::matrixTimesMatrixELL(
    const ValueType alpha,
    const ELLStorage<ValueType>& a,
    const ELLStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B, with " << "A = " << a << ", B = " << b << ", all are ELL" )

    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<ELLKernelTrait::matrixMultiplySizes> matrixMultiplySizes;
    static LAMAKernel<ELLKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    ContextPtr loc = Context::getHostPtr();  // not yet available on other devices

    SCAI_ASSERT_ERROR( &a != this, "matrixTimesMatrix: alias of a with this result matrix" )
    SCAI_ASSERT_ERROR( &b != this, "matrixTimesMatrix: alias of b with this result matrix" )

    SCAI_ASSERT_EQUAL_ERROR( a.getNumColumns(), b.getNumRows() )

    allocate( a.getNumRows(), b.getNumColumns() );

    mDiagonalProperty = ( mNumRows == mNumColumns );

    {
        ReadAccess<IndexType> aIA( a.getIA(), loc );
        ReadAccess<IndexType> aJA( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );

        ReadAccess<IndexType> bIA( b.getIA(), loc );
        ReadAccess<IndexType> bJA( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );

        allocate( a.getNumRows(), b.getNumColumns() );

        WriteOnlyAccess<IndexType> cIA( mIA, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        // 1. Step: compute resulting IA array
        matrixMultiplySizes[loc] ( cIA.get(), a.getNumRows(), a.getNumColumns(), b.getNumRows(), false, aIA.get(), aJA.get(),
                                   a.getNumValuesPerRow(), bIA.get(), bJA.get(), b.getNumValuesPerRow() );

        // 2. Step: compute length of longest row
        mNumValuesPerRow = reduce[ loc ]( cIA.get(), mNumRows, common::reduction::MAX );

        // 3. Step: Allocate IA and Values arrays with new size
        WriteOnlyAccess<IndexType> cJA( mJA, loc, mNumValuesPerRow * mNumRows );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValuesPerRow * mNumRows );

        // 4. Step: Compute cJA and cValues
        matrixMultiply[loc]( cJA.get(), cValues.get(), cIA.get(), mNumValuesPerRow, mNumRows, mNumColumns, b.getNumRows(),
                             false, alpha, aIA.get(), aJA.get(), aValues.get(), a.getNumValuesPerRow(), bIA.get(), bJA.get(),
                             bValues.get(), b.getNumValuesPerRow() );
    }

    // 5. Step: Computation of C might have produced some zero elements

    compress();
}

template<typename ValueType>
void ELLStorage<ValueType>::matrixAddMatrixELL(
    const ValueType alpha,
    const ELLStorage<ValueType>& a,
    const ValueType beta,
    const ELLStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " * A + " << beta << " * B, with " << "A = " << a << ", B = " << b << ", all are ELL" )

    static LAMAKernel<ELLKernelTrait::matrixAddSizes> matrixAddSizes;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<ELLKernelTrait::matrixAdd<ValueType> > matrixAdd;

    ContextPtr loc = matrixAddSizes.getValidContext( reduce, matrixAdd, this->getContextPtr() );

    SCAI_ASSERT_ERROR( &a != this, "matrixAddMatrix: alias of a with this result matrix" )
    SCAI_ASSERT_ERROR( &b != this, "matrixAddMatrix: alias of b with this result matrix" )

    allocate( a.getNumRows(), a.getNumColumns() );

    SCAI_ASSERT_EQUAL_ERROR( mNumRows, b.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, b.getNumColumns() )

    //mDiagonalProperty = ( mNumRows == mNumColumns );

    {
        ReadAccess<IndexType> aIA( a.getIA(), loc );
        ReadAccess<IndexType> aJA( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );

        ReadAccess<IndexType> bIA( b.getIA(), loc );
        ReadAccess<IndexType> bJA( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );

        WriteOnlyAccess<IndexType> cIA( mIA, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        // 1. Step: Compute IA array
        matrixAddSizes[loc]( cIA.get(), a.getNumRows(), a.getNumColumns(), false, aIA.get(), aJA.get(),
                             a.getNumValuesPerRow(), bIA.get(), bJA.get(), b.getNumValuesPerRow() );

        // 2. Step: compute length of longest row
        mNumValuesPerRow = reduce[loc]( cIA.get(), mNumRows, common::reduction::MAX );

        // 3. Step: Allocate IA and Values arrays with new size
        WriteOnlyAccess<IndexType> cJA( mJA, loc, mNumValuesPerRow * mNumRows );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValuesPerRow * mNumRows );

        // 4. Step: Compute cJA and cValues
        matrixAdd[loc]( cJA.get(), cValues.get(), cIA.get(), mNumValuesPerRow, mNumRows, mNumColumns, false, alpha,
                        aIA.get(), aJA.get(), aValues.get(), a.getNumValuesPerRow(), beta, bIA.get(), bJA.get(),
                        bValues.get(), b.getNumValuesPerRow() );

    }

    // 5. Step: Computation of C might have produced some zero elements
    compress();

    check( "result of matrix + matrix" ); // just verify for a correct matrix
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>* ELLStorage<ValueType>::copy() const
{
    SCAI_LOG_INFO( logger, "copy" )

    return new ELLStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>* ELLStorage<ValueType>::newMatrixStorage() const
{
   common::unique_ptr<ELLStorage<ValueType> > storage( new ELLStorage<ValueType>() ); 
   storage->setContextPtr( this->getContextPtr() );
   return storage.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* ELLStorage<ValueType>::create()
{
    return new ELLStorage<ValueType>();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType ELLStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::ELL, common::getScalarType<ValueType>() );
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_ELL_STORAGE_INSTANTIATE(z, I, _)                                     \
                                                                                  \
    template<>                                                                    \
    const char* ELLStorage<ARITHMETIC_HOST_TYPE_##I>::typeName()                  \
    {                                                                             \
        return "ELLStorage<" PRINT_STRING(ARITHMETIC_HOST_TYPE_##I) ">";      \
    }                                                                             \
                                                                                  \
    template class COMMON_DLL_IMPORTEXPORT ELLStorage<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_ELL_STORAGE_INSTANTIATE, _ )

#undef LAMA_ELL_STORAGE_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
