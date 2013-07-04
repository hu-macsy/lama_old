/**
 * @file ELLStorage.cpp
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
 * @brief Instantitions for template class ELLStorage.
 * @author Lauretta Schubert
 * @date 25.05.2011
 * @since 1.0.0
 */

// hpp
#include <lama/storage/ELLStorage.hpp>

// others
#include <lama/ContextAccess.hpp>
#include <lama/HostReadAccess.hpp>
#include <lama/LAMAArrayUtils.hpp>
#include <lama/LAMAInterface.hpp>
#include <lama/task/TaskSyncToken.hpp>
#include <lama/NoSyncToken.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/WriteAccess.hpp>

// tracing
#include <lama/tracing.hpp>

namespace lama
{

using boost::shared_ptr;

/* --------------------------------------------------------------------------- */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ELLStorage<ValueType>::logger, "MatrixStorage.ELLStorage" )

/* --------------------------------------------------------------------------- */

template<>
const char* ELLStorage<float>::typeName()
{
    return "ELLStorage<float>";
}

template<>
const char* ELLStorage<double>::typeName()
{
    return "ELLStorage<double>";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const ContextType con /* = Context::Host */)

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>( numRows, numColumns ), mNumValuesPerRow( 0 )
{
    // TODO in other formats the last parameter is "const ContextPtr loc"

    ContextPtr loc = ContextFactory::getContext( con );
    setContext( loc );

    // Initialization requires correct values for the IA array with 0

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )

    WriteOnlyAccess<IndexType> ellSizes( mIA, loc, mNumRows );

    LAMA_CONTEXT_ACCESS( loc )

    setVal( ellSizes.get(), mNumRows, 0 );

    LAMA_LOG_DEBUG( logger, "ELLStorage for matrix " << mNumRows << " x " << mNumColumns << ", no elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage()

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>( 0, 0 ), mNumValuesPerRow( 0 )
{
    LAMA_LOG_DEBUG( logger, "ELLStorage, default constructor for zero matrix." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRows,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<ValueType>& values )

    : CRTPMatrixStorage<ELLStorage<ValueType>,ValueType>()
{
    LAMA_LOG_INFO( logger, "constructor with ELL data array" )

    setELLData( numRows, numColumns, numValuesPerRows, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( const ELLStorage<ValueType>& other )

    : CRTPMatrixStorage<ELLStorage<ValueType>, ValueType>( 0, 0 )
{
    LAMA_LOG_INFO( logger, "constructor # other = " << other )

    // call the assignment operator

    operator=( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>& ELLStorage<ValueType>::operator=( const ELLStorage<ValueType>& other )
{
    // nothing to do for self assignments

    if ( &other == this )
    {
        return *this;
    }

    // For optimization: fillValues, resetDiagonalProperty, check is not really needed

    setELLData( other.mNumRows, other.mNumColumns, other.mNumValuesPerRow,
                other.mIA, other.mJA, other.mValues );

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>& ELLStorage<ValueType>::operator=( const _MatrixStorage& other )
{
    assign( other );   // calls virtual method of MatrixStorage

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::print() const
{
    LAMA_LOG_INFO( logger, "print" )

    using std::cout;
    using std::endl;

    cout << "ELLStorage " << mNumRows << " x " << mNumColumns << ", #values = " << getNumValues() << endl;

    HostReadAccess<IndexType> ia( mIA );
    HostReadAccess<IndexType> ja( mJA );
    HostReadAccess<ValueType> values( mValues );

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        cout << "Row " << i << " ( " << ia[i] << " entries ) :";

        for ( IndexType jj = 0; jj < ia[i]; ++jj )
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
    LAMA_LOG_INFO( logger, "getFormat" )
    return ELL;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType ELLStorage<ValueType>::getNumValues() const
{
    LAMA_LOG_INFO( logger, "getNumValues" )

    const ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( sum, loc, Utils, Reductions, IndexType )

    ReadAccess<IndexType> ia( mIA, loc );

    LAMA_CONTEXT_ACCESS( loc )

    IndexType numValues = sum( ia.get(), mNumRows );

    return numValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::purge()
{
    LAMA_LOG_INFO( logger, "purge" )

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
    LAMA_LOG_INFO( logger, "set identity # size = " << size )

    ContextPtr loc = getContextPtr();

    mNumRows = size;
    mNumColumns = size;
    mNumValuesPerRow = 1;

    {
        WriteOnlyAccess<IndexType> ia( mIA, loc, mNumRows );
        WriteOnlyAccess<IndexType> ja( mJA, loc, mNumRows );

        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )
        LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )

        LAMA_CONTEXT_ACCESS( loc )

        setVal( ia.get(), mNumRows, 1 );
        setOrder( ja.get(), mNumRows );
    }

    // extra block caused by differnt types of setVal()
    {
        WriteOnlyAccess<ValueType> data( mValues, loc, mNumRows );

        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

        LAMA_CONTEXT_ACCESS( loc )

        ValueType one = static_cast<ValueType>( 1.0 );
        setVal( data.get(), mNumRows, one );
    }

    mDiagonalProperty = true;

    LAMA_LOG_INFO( logger, *this << " is identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool ELLStorage<ValueType>::checkDiagonalProperty() const
{
    LAMA_LOG_INFO( logger, "checkDiagonalProperty" )

    IndexType numDiagonals = std::min( mNumRows, mNumColumns );

    bool diagonalProperty = true;

    if ( numDiagonals == 0 )
    {
        // diagonal property is given for zero-sized matrices

        diagonalProperty = true;
    }
    else if ( mNumValuesPerRow < 1 )
    {
        // no elements, so certainly it does not have diagonl property

        diagonalProperty = false;
    }
    else
    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN( hasDiagonalProperty, loc, ELLUtils, Operations )

        ReadAccess<IndexType> ja( mJA, loc );

        LAMA_CONTEXT_ACCESS( loc )

        diagonalProperty = hasDiagonalProperty( numDiagonals, ja.get() );
    }

    LAMA_LOG_INFO( logger, *this << ": checkDiagonalProperty = " << diagonalProperty )

    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::clear()
{
    LAMA_LOG_INFO( logger, "clear" )

    mNumRows         = 0;
    mNumColumns      = 0;
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
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr loc ) const
{
    LAMA_REGION( "Storage.ELL->CSR" )

    LAMA_LOG_INFO( logger,
                   "buildCSR<" << Scalar::getType<OtherValueType>() << ">" << " from ELL<" << Scalar::getType<ValueType>() << ">" << " on " << *loc )

    LAMA_INTERFACE_FN( sizes2offsets, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, IndexType, IndexType )
    LAMA_INTERFACE_FN_TT( getCSRValues, loc, ELLUtils, Conversions, ValueType, OtherValueType )

    ReadAccess<IndexType> ellSizes( mIA, loc );
    WriteAccess<IndexType> csrIA( ia, loc );

    csrIA.resize( mNumRows + 1 );

    LAMA_CONTEXT_ACCESS( loc )
    // just copy the size array mIA
    set( csrIA.get(), ellSizes.get(), mNumRows );

    if ( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets( csrIA.get(), mNumRows );

    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );

    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<OtherValueType> csrValues( *values, loc, numValues );

    getCSRValues( csrJA.get(), csrValues.get(), csrIA.get(), mNumRows, mNumValuesPerRow, ellSizes.get(), ellJA.get(), ellValues.get() );

}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    const ContextPtr loc )
{
    LAMA_REGION( "Storage.ELL<-CSR" )

    LAMA_LOG_INFO( logger, "set CSR data on " << *loc 
                            << ": numRows = " << numRows << ", numColumns = " << numColumns 
                            << ", numValues = " << numValues 
                            << ", compress threshold = " << mCompressThreshold )

    if ( numRows == 0 )
    {
        // just allocate will clear member arrays

        allocate( numRows, numColumns );

        return;
    }

    _MatrixStorage::setDimension( numRows, numColumns );

    // Get function pointers for needed routines at the LAMA interface

    LAMA_INTERFACE_FN( offsets2sizes, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN( hasDiagonalProperty, loc, ELLUtils, Operations )
    LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )
    LAMA_INTERFACE_FN_TT( setCSRValues, loc, ELLUtils, Conversions, ValueType, OtherValueType )

    // build array with non-zero values per row

    {
        ReadAccess<IndexType> csrIA( ia, loc );
        WriteOnlyAccess<IndexType> ellSizes( mIA, loc, mNumRows );

        LAMA_CONTEXT_ACCESS( loc )
        offsets2sizes( ellSizes.get(), csrIA.get(), mNumRows );
    }

    // determine the maximal number of non-zero in one row

    {
        ReadAccess<IndexType> ellSizes( mIA, loc );
        LAMA_CONTEXT_ACCESS( loc )
        mNumValuesPerRow = maxval( ellSizes.get(), mNumRows );
    }

    LAMA_LOG_INFO( logger, "setCSRData, #values/row = " << mNumValuesPerRow )

    //  Now we know the size of the ja and values arrays for the ELL format

    const IndexType dataSize = mNumValuesPerRow * mNumRows;

    if ( mNumRows > 200 && mNumValuesPerRow > 0 )
    {
        // make this check only on larger matrices, dataSize must not be equal 0

        double fillRate = double( numValues ) / double( dataSize );

        if ( fillRate < 0.5 )
        {
            LAMA_LOG_WARN( logger, *this << ": fill rate = " << fillRate 
                           << " ( " << numValues << " non-zero values ), consider using JDS" )
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

        LAMA_LOG_DEBUG( logger, "convert CSR -> ELL, ellSize = " << dataSize )

        LAMA_CONTEXT_ACCESS( loc )

        setCSRValues( ellJA.get(), ellValues.get(), 
                      ellIA.get(), mNumRows, mNumValuesPerRow, 
                      csrIA.get(), csrJA.get(), csrValues.get() 
                    );

        LAMA_LOG_DEBUG( logger, " size = " <<ellJA.size() )

        IndexType numDiagonals = std::min( mNumRows, mNumColumns );

        if ( numDiagonals == 0)
        { 
            mDiagonalProperty = true;
        }
        else if ( numValues == 0 )
        {
            mDiagonalProperty = false;
        } 
        else
        {
            LAMA_CONTEXT_ACCESS( loc )
            mDiagonalProperty = hasDiagonalProperty( numDiagonals, ellJA.get() );
        }
    }

    if ( numRows == numColumns && !mDiagonalProperty )
    {
        LAMA_LOG_INFO( logger, *this << ": square matrix has not diagonal property" )
    }

    buildRowIndexes( loc );

    LAMA_LOG_DEBUG( logger, "convert CSR -> ELL done: " << *this )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void ELLStorage<ValueType>::setELLData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const _LAMAArray& values )
{
    LAMA_ASSERT_EQUAL_ERROR( numRows, ia.size() )
    LAMA_ASSERT_EQUAL_ERROR( numRows * numValuesPerRow, ja.size() )
    LAMA_ASSERT_EQUAL_ERROR( numRows * numValuesPerRow, values.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValuesPerRow = numValuesPerRow;

    ContextPtr loc = getContextPtr();

    LAMAArrayUtils::assignImpl( mIA, ia, loc );
    LAMAArrayUtils::assignImpl( mJA, ja, loc );

    LAMAArrayUtils::assign( mValues, values, loc );  // supports type conversion

    // fill up my arrays ja and values to make matrix-multiplication fast

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT_T( fillELLValues, loc, ELLUtils, Solver, ValueType )

        ReadAccess<IndexType> ellIA( mIA, loc );
        WriteAccess<IndexType> ellJA( mJA, loc );
        WriteAccess<ValueType> ellValues( mValues, loc );

        LAMA_LOG_DEBUG( logger, "fill ELL data" )

        LAMA_CONTEXT_ACCESS( loc )

        fillELLValues( ellJA.get(), ellValues.get(),
                       ellIA.get(), mNumRows, mNumValuesPerRow );
    }

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode

#ifdef LAMA_ASSERT_LEVEL_DEBUG
    check( "ELLStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif

    this->resetDiagonalProperty();

    LAMA_LOG_INFO( logger, *this << ": set ELLPACK by arrays ia, ja, values" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType ELLStorage<ValueType>::getNumValuesPerRow() const
{
    LAMA_LOG_INFO( logger, "getNumValuesPerRow" )
    return mNumValuesPerRow;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setDiagonalImpl( const Scalar scalar )
{
    LAMA_LOG_INFO( logger, "setDiagonalImpl # scalar = " << scalar )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    ValueType value = scalar.getValue<ValueType>();

    LAMA_CONTEXT_ACCESS( loc )
    setVal( wValues.get(), numDiagonalElements, value );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherType>& diagonal )
{
    LAMA_LOG_INFO( logger, "setDiagonalImpl # diagonal = " << diagonal )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, ValueType, OtherType )

    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    // ELL format with diagonal property: diagonal is just the first column in mValues

    LAMA_CONTEXT_ACCESS( loc )

    set( wValues.get(), rDiagonal.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::getRowImpl( LAMAArray<OtherType>& row, const IndexType i ) const
{
    LAMA_LOG_TRACE( logger, "getRowImpl # row = " << row << ", i = " << i )

    LAMA_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( getRow, loc, ELLUtils, Getter, ValueType, OtherType )

    WriteOnlyAccess<OtherType> wRow( row, loc, mNumColumns );
    const ReadAccess<IndexType> rIa( mIA, loc );
    const ReadAccess<IndexType> rJa( mJA, loc );
    const ReadAccess<ValueType> rValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    getRow( wRow.get(), i, mNumRows, mNumColumns, mNumValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherType>& diagonal ) const
{
    LAMA_LOG_INFO( logger, "getDiagonalImpl # diagonal = " << diagonal )

    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( set, loc, Utils, Copy, OtherType, ValueType )

    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );

    // ELL format with diagonal property: diagonal is just the first column in mValues

    LAMA_CONTEXT_ACCESS( loc )

    set( wDiagonal.get(), rValues.get(), numDiagonalElements );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    LAMA_LOG_INFO( logger, "scaleImpl # scalar = " << scalar )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( scale, loc, Utils, Transform, ValueType )

    WriteAccess<ValueType> wValues( mValues, loc );

    const ValueType value = scalar.getValue<ValueType>();

    LAMA_CONTEXT_ACCESS( loc )

    scale( wValues.get(), value, mValues.size() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::scaleImpl( const LAMAArray<OtherValueType>& values )
{
    LAMA_LOG_INFO( logger, "scaleImpl # values = " << values )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( scaleValue, loc, ELLUtils, Scale, ValueType, OtherValueType )

    ReadAccess<OtherValueType> rValues( values, loc );
    ReadAccess<IndexType> rIa( mIA, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    scaleValue( mNumRows, mNumValuesPerRow, rIa.get(), wValues.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& ELLStorage<ValueType>::getIA() const
{
    return mIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<IndexType>& ELLStorage<ValueType>::getJA() const
{
    return mJA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const LAMAArray<ValueType>& ELLStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::~ELLStorage()
{
    LAMA_LOG_DEBUG( logger,
                    "~ELLStorage for matrix " << mNumRows << " x " << mNumColumns << ", # nnr = " << mNumValuesPerRow )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::check( const char* msg ) const
{
    LAMA_LOG_INFO( logger, "check # msg = " << msg )

    LAMA_ASSERT_EQUAL_ERROR( mNumRows, mIA.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValuesPerRow * mNumRows, mJA.size() )
    LAMA_ASSERT_EQUAL_ERROR( mJA.size(), mValues.size() )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN( check, loc, ELLUtils, Operations )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    LAMA_CONTEXT_ACCESS( loc )

    check( mNumRows, mNumValuesPerRow, mNumColumns, rIa.get(), rJa.get(), msg );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    LAMA_LOG_INFO( logger, "allocate ELL sparse matrix of size " << numRows << " x " << numColumns )

    clear();   

    mNumRows    = numRows;
    mNumColumns = numColumns;

    LAMA_LOG_DEBUG( logger, "resize mIA, mNumRows = " << mNumRows )

    {
        // Intialize array mIA with 0

        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )

        WriteOnlyAccess<IndexType> ia( mIA, loc, mNumRows );

        LAMA_CONTEXT_ACCESS( loc )

        setVal( ia.get(), mNumRows, 0 );
    }

    mDiagonalProperty = checkDiagonalProperty();

    LAMA_LOG_DEBUG( logger, "ready allocate" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ELLStorage<" << Scalar::getType<ValueType>() << ">( rows=" << mNumRows << ", cols=" << mNumColumns
           << ", nnr=" << mNumValuesPerRow << ", threshold = " << mCompressThreshold << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    LAMA_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )
    LAMA_LOG_TRACE( logger, "sizes: ia = " << mIA.size() << ", ja = " << mJA.size() << ", data = " << mValues.size() )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( getValue, loc, ELLUtils, Getter, ValueType, ValueType )

    const ReadAccess<IndexType> rIa( mIA, loc );
    const ReadAccess<IndexType> rJa( mJA, loc );
    const ReadAccess<ValueType> rValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    return getValue( i, j, mNumRows, mNumValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    LAMA_LOG_INFO( logger, "prefetch # location " << location )
    LAMA_LOG_DEBUG( logger, "Starting prefetch of "<<*this<<" to "<<location )

    mRowIndexes.prefetch( location );
    mIA.prefetch( location );
    mJA.prefetch( location );
    mValues.prefetch( location );

    LAMA_LOG_DEBUG( logger, "Finished prefetch of "<<*this<<" to "<<location )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::wait() const
{
    LAMA_LOG_INFO( logger, "wait" )

    mRowIndexes.wait();
    mIA.wait();
    mJA.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::buildRowIndexes( const ContextPtr loc )
{
    LAMA_LOG_INFO( logger, "buildRowIndexes # loc = " << loc )

    mRowIndexes.clear();

    if ( mNumRows == 0 )
    {
        return;
    }

    // Get function pointers for needed routines at the LAMA interface
    LAMA_INTERFACE_FN( countNonEmptyRowsBySizes, loc, ELLUtils, Operations )
    LAMA_INTERFACE_FN( setNonEmptyRowsBySizes, loc, ELLUtils, Operations )

    ReadAccess<IndexType> ellIA( mIA, loc );
    LAMA_CONTEXT_ACCESS( loc )
    IndexType nonZeroRows = countNonEmptyRowsBySizes( ellIA.get(), mNumRows );

    float usage = float( nonZeroRows ) / float( mNumRows );

    if ( usage >= mCompressThreshold )
    {
        LAMA_LOG_INFO( logger, "ELLStorage: do not build row indexes, usage = " << usage 
                                << " >= " << mCompressThreshold << " ( threshold )" )
        return;
    }

    WriteOnlyAccess<IndexType> rowIndexes( mRowIndexes, loc, nonZeroRows );
    setNonEmptyRowsBySizes( rowIndexes.get(), nonZeroRows, ellIA.get(), mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::compress( const ValueType eps /* = 0.0 */)
{
    LAMA_LOG_INFO( logger, "compress: eps = " << eps )

    // TODO: Implement for CUDA
    ContextPtr loc = ContextFactory::getContext( Context::Host );

    LAMA_INTERFACE_FN_T( compressIA, loc, ELLUtils, Helper, ValueType )
    LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )
    LAMA_INTERFACE_FN_T( compressValues, loc, ELLUtils, Helper, ValueType )

    ReadAccess<IndexType> IA( mIA, loc );
    ReadAccess<IndexType> JA( mJA, loc );
    ReadAccess<ValueType> values( mValues, loc );

    // 1. Step: Check for 0 elements and write new IA array
    LAMAArray<IndexType> newIAArray;
    WriteOnlyAccess<IndexType> newIA( newIAArray, loc, mNumRows );

    compressIA( IA.get(), JA.get(), values.get(), mNumRows, mNumValuesPerRow, eps, newIA.get() );

    // 2. Step: compute length of longest row
    IndexType newNumValuesPerRow = maxval( IA.get(), mNumRows );

    // Do further steps, if new array could be smaller
    if ( newNumValuesPerRow < mNumValuesPerRow )
    {
        // 3. Step: Allocate new JA and Values array
        LAMAArray<ValueType> newValuesArray;
        LAMAArray<IndexType> newJAArray;
        WriteOnlyAccess<ValueType> newValues( newValuesArray, loc, mNumRows * newNumValuesPerRow );
        WriteOnlyAccess<IndexType> newJA( newJAArray, loc, mNumRows * newNumValuesPerRow );

        // 4. Step: Compute new JA and Values array
        compressValues( IA.get(), JA.get(), values.get(), mNumRows, mNumValuesPerRow, eps, newNumValuesPerRow, newJA.get(), newValues.get() );

        mJA.swap( newJAArray );
        mValues.swap( newValuesArray );
        mNumValuesPerRow = newNumValuesPerRow;
    }
}

template<typename ValueType>
void ELLStorage<ValueType>::swap( ELLStorage<ValueType>& other )
{
    LAMA_LOG_INFO( logger, "swap # other = " << other )

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
    LAMA_LOG_INFO( logger, "getMemoryUsageImpl" )

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
    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const
{
    LAMA_LOG_INFO( logger,
                   *this << ": matrixTimesVector, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    ContextPtr loc = getContextPtr();

    if ( mNumValuesPerRow == 0 )
    {
        // this matrix is ZERO, so all to do is result = beta * y

        LAMAArrayUtils::assignScaled( result, beta, y, loc );
        return;
    }

    LAMA_REGION( "Storage.ELL.timesVector" )

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_T( sparseGEMV, loc, ELLUtils, Mult, ValueType )
    LAMA_INTERFACE_FN_T( normalGEMV, loc, ELLUtils, Mult, ValueType )

    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );

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
            sparseGEMV( wResult.get(), mNumRows, mNumValuesPerRow, alpha, rX.get(), numNonZeroRows, rows.get(),
                        ellIA.get(), ellJA.get(), ellValues.get(), NULL );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            LAMA_CONTEXT_ACCESS( loc )
            normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, mNumValuesPerRow, ellIA.get(),
                        ellJA.get(), ellValues.get(), NULL );
        }
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc )
        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumValuesPerRow, ellIA.get(),
                    ellJA.get(), ellValues.get(), NULL );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::matrixTimesVectorAsync(
    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const
{
    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    LAMA_REGION( "Storage.ELL.timesVectorAsync" )

    ContextPtr loc = getContextPtr();

    if ( mNumValuesPerRow == 0 )
    {
        // this matrix is ZERO, so all to do is result = beta * y
        // that is already sychronized here

        LAMAArrayUtils::assignScaled( result, beta, y, loc );
        return new NoSyncToken();
    }

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVectorAsync on " << *loc )

    LAMA_INTERFACE_FN_T( sparseGEMV, loc, ELLUtils, Mult, ValueType )
    LAMA_INTERFACE_FN_T( normalGEMV, loc, ELLUtils, Mult, ValueType )

    std::auto_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > ellIA( new ReadAccess<IndexType>( mIA, loc ) );
    shared_ptr<ReadAccess<IndexType> > ellJA( new ReadAccess<IndexType>( mJA, loc ) );
    shared_ptr<ReadAccess<ValueType> > ellValues( new ReadAccess<ValueType>( mValues, loc ) );
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

            shared_ptr<ReadAccess<IndexType> > rRowIndexes( new ReadAccess<IndexType>( mRowIndexes, loc ) );

            syncToken->pushAccess( rRowIndexes );

            LAMA_CONTEXT_ACCESS( loc )

            sparseGEMV( wResult->get(), mNumRows, mNumValuesPerRow, alpha, rX->get(), numNonZeroRows, rRowIndexes->get(),
                        ellIA->get(), ellJA->get(), ellValues->get(), syncToken.get() );
        }
        else
        {
            // we assume that normalGEMV can deal with the alias of result, y

            LAMA_CONTEXT_ACCESS( loc )

            normalGEMV( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, mNumValuesPerRow,
                        ellIA->get(), ellJA->get(), ellValues->get(), syncToken.get() );
        }

        syncToken->pushAccess( wResult );
    }
    else
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumRows ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        LAMA_CONTEXT_ACCESS( loc )

        normalGEMV( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, mNumValuesPerRow, ellIA->get(),
                    ellJA->get(), ellValues->get(), syncToken.get() );

        syncToken->pushAccess( wResult );
        syncToken->pushAccess( rY );
    }

    syncToken->pushAccess( ellIA );
    syncToken->pushAccess( ellJA );
    syncToken->pushAccess( ellValues );
    syncToken->pushAccess( rX );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.ELL.jacobiIterate" )

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

    LAMA_INTERFACE_FN_T( jacobi, loc, ELLUtils, Solver, ValueType )

    // make all needed data available at loc

    WriteAccess<ValueType> wSolution( solution, loc );
    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );

    LAMA_CONTEXT_ACCESS( loc )

    jacobi( wSolution.get(), mNumRows, mNumValuesPerRow, ellSizes.get(), ellJA.get(), ellValues.get(),
            rOldSolution.get(), rRhs.get(), omega, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::jacobiIterateAsync(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.ELL.jacobiIterateAsync" )

    ContextPtr loc = getContextPtr();

    if ( loc->getType() == Context::Host )
    {
        // used later in OpenMP to generate a TaskSyncToken

        void ( ELLStorage::*jb )(
            LAMAArray<ValueType>&,
            const LAMAArray<ValueType>&,
            const LAMAArray<ValueType>&,
            const ValueType omega ) const

        = &ELLStorage<ValueType>::jacobiIterate;

        using boost::bind;
        using boost::cref;
        using boost::ref;

        return new TaskSyncToken( bind( jb, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task

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

    LAMA_INTERFACE_FN_T( jacobi, loc, ELLUtils, Solver, ValueType )

    std::auto_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // make all needed data available at loc

    shared_ptr<WriteAccess<ValueType> > wSolution( new WriteAccess<ValueType>( solution, loc ) );
    shared_ptr<ReadAccess<IndexType> > ellSizes( new ReadAccess<IndexType>( mIA, loc ) );
    shared_ptr<ReadAccess<IndexType> > ellJA( new ReadAccess<IndexType>( mJA, loc ) );
    shared_ptr<ReadAccess<ValueType> > ellValues( new ReadAccess<ValueType>( mValues, loc ) );
    shared_ptr<ReadAccess<ValueType> > rOldSolution( new ReadAccess<ValueType>( oldSolution, loc ) );
    shared_ptr<ReadAccess<ValueType> > rRhs( new ReadAccess<ValueType>( rhs, loc ) );

    LAMA_CONTEXT_ACCESS( loc )

    jacobi( wSolution->get(), mNumRows, mNumValuesPerRow, ellSizes->get(), ellJA->get(), ellValues->get(),
            rOldSolution->get(), rRhs->get(), omega, syncToken.get() );

    syncToken->pushAccess( rRhs );
    syncToken->pushAccess( rOldSolution );
    syncToken->pushAccess( ellValues );
    syncToken->pushAccess( ellJA );
    syncToken->pushAccess( ellSizes );
    syncToken->pushAccess( wSolution );

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.ELL.jacobiIterateHalo" )

    LAMA_LOG_INFO( logger, "HOST: Jacobi iteration on halo matrix data." )

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumRows() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumColumns() )
    LAMA_ASSERT_DEBUG( localStorage.hasDiagonalProperty(), localStorage << ": has not diagonal property" )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    const LAMAArray<ValueType>* localDiagonal;

    // might be we need a temporary LAMA array for the local diagonal

    boost::shared_ptr<LAMAArray<ValueType> > tmpLocalDiagonal;

    if ( localStorage.getFormat() == ELL )
    {
        const ELLStorage<ValueType>* ellLocal;

        ellLocal = dynamic_cast<const ELLStorage<ValueType>*>( &localStorage );
        LAMA_ASSERT_DEBUG( ellLocal, "could not cast to ELLStorage " << localStorage )
        localDiagonal = &( ellLocal->mValues );
    }
    else
    {
        // make a temporary for the diagonal and get it from local storage

        LAMA_LOG_WARN( logger, "local stroage is not ELL, temorary needed for diagonal" )

        tmpLocalDiagonal = boost::shared_ptr<LAMAArray<ValueType> >( new LAMAArray<ValueType>() );
        localStorage.getDiagonal( *tmpLocalDiagonal );
        localDiagonal = tmpLocalDiagonal.get();

        // Note: tmpLocalDiagonal will be freed at end of routine
    }

    jacobiIterateHalo( localSolution, localDiagonal, oldHaloSolution, omega );

}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const LAMAArray<ValueType>* localDiagonal,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.ELL.jacobiIterateHalo" )

    LAMA_LOG_INFO( logger, "HOST: Jacobi iteration on halo matrix data." )

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobiHalo, loc, ELLUtils, Solver, ValueType )

    {
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<ValueType> rLocalDiagonal( *localDiagonal, loc );
        ReadAccess<IndexType> haloIA( mIA, loc );
        ReadAccess<IndexType> haloJA( mJA, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

        const IndexType numNonEmptyRows = mRowIndexes.size();

        if ( numNonEmptyRows != 0 )
        {
            ReadAccess<IndexType> haloRowIndexes( mRowIndexes, loc );

            LAMA_CONTEXT_ACCESS( loc )

            jacobiHalo( wSolution.get(), mNumRows, rLocalDiagonal.get(), mNumValuesPerRow, haloIA.get(), haloJA.get(),
                        haloValues.get(), haloRowIndexes.get(), numNonEmptyRows, rOldHaloSolution.get(), omega, NULL );
        }
        else
        {
            // no row indexes available, computation is done over all rows

            const IndexType numNonEmptyRows = mNumRows;

            LAMA_CONTEXT_ACCESS( loc )

            jacobiHalo( wSolution.get(), mNumRows, rLocalDiagonal.get(), mNumValuesPerRow, haloIA.get(), haloJA.get(),
                        haloValues.get(), NULL, numNonEmptyRows, rOldHaloSolution.get(), omega, NULL );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::maxNorm() const
{
    LAMA_LOG_INFO( logger, *this << ": maxNorm()" )

    if ( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return 0.0f;
    }

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, ELLUtils, Reductions, ValueType )

    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal( mNumRows, mNumValuesPerRow, ellIA.get(), ellValues.get() );

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
    LAMA_LOG_INFO( logger,
                   "this = " << alpha << " * A * B + " << beta << " * C, with " << "A = " << a << ", B = " << b << ", C = " << c )

    const ELLStorage<ValueType>* ellA = NULL;
    const ELLStorage<ValueType>* ellB = NULL;
    const ELLStorage<ValueType>* ellC = NULL;

    //    boost::shared_ptr<CSRStorage<ValueType> > tmpA;
    //    boost::shared_ptr<CSRStorage<ValueType> > tmpB;
    boost::shared_ptr<ELLStorage<ValueType> > tmpC;

    if ( a.getFormat() == ELL )
    {
        ellA = dynamic_cast<const ELLStorage<ValueType>*>( &a );
        LAMA_ASSERT_DEBUG( ellA, "could not cast to ELLStorage " << a )
    }
    else
    {
        LAMA_LOG_ERROR( logger, a << ": a not ELL format" )
    }
    if ( b.getFormat() == ELL )
    {
        ellB = dynamic_cast<const ELLStorage<ValueType>*>( &b );
        LAMA_ASSERT_DEBUG( ellB, "could not cast to ELLStorage " << b )
    }
    else
    {
        LAMA_UNSUPPORTED( b << ": b not ELL format" )
    }

    if ( ellA == NULL || ellB == NULL )
    {
        // input matrices not ELL format, so try via CSR

        MatrixStorage<ValueType>::matrixTimesMatrix( alpha, a, b, beta, c );
        return;
    }

    if ( beta != 0.0 )
    {
        if ( ( c.getFormat() == ELL ) && ( &c != this ) )
        {
            ellC = dynamic_cast<const ELLStorage<ValueType>*>( &c );
            LAMA_ASSERT_DEBUG( ellC, "could not cast to ELLStorage " << c )
        }
        else
        {
            LAMA_UNSUPPORTED( c << ": ELL temporary required for matrix add" )
            tmpC = boost::shared_ptr<ELLStorage<ValueType> >( new ELLStorage<ValueType>( c ) );
            ellC = tmpC.get();
        }

    }

    ELLStorage<ValueType> tmp;
    tmp.matrixTimesMatrixELL( alpha, *ellA, *ellB );

    if ( beta != 0 )
    {
        ELLStorage<ValueType> tmp1;
        tmp1.matrixAddMatrixELL( 1.0, tmp, beta, *ellC );
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
    LAMA_LOG_INFO( logger,
                   *this << ": = " << alpha << " * A * B, with " << "A = " << a << ", B = " << b << ", all are ELL" )

    // TODO: Implement for CUDA
    ContextPtr loc = ContextFactory::getContext( Context::Host );

    LAMA_INTERFACE_FN( matrixMultiplySizes, loc, ELLUtils, MatrixExpBuild )
    LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )
    LAMA_INTERFACE_FN_T( matrixMultiply, loc, ELLUtils, MatrixExp, ValueType )

    LAMA_ASSERT_ERROR( &a != this, "matrixTimesMatrix: alias of a with this result matrix" )
    LAMA_ASSERT_ERROR( &b != this, "matrixTimesMatrix: alias of b with this result matrix" )

    LAMA_ASSERT_EQUAL_ERROR( a.getNumColumns(), b.getNumRows() )

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

        LAMA_CONTEXT_ACCESS( loc )

        // 1. Step: compute resulting IA array
        matrixMultiplySizes( cIA.get(),
                             a.getNumRows(), a.getNumColumns(), b.getNumRows(), false,
                             aIA.get(), aJA.get(), a.getNumValuesPerRow(),
                             bIA.get(), bJA.get(), b.getNumValuesPerRow() );

        // 2. Step: compute length of longest row
        mNumValuesPerRow = maxval( cIA.get(), mNumRows );

        // 3. Step: Allocate IA and Values arrays with new size
        WriteOnlyAccess<IndexType> cJA( mJA, loc, mNumValuesPerRow * mNumRows );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValuesPerRow * mNumRows );

        // 4. Step: Compute cJA and cValues
        matrixMultiply( cJA.get(), cValues.get(), cIA.get(), mNumValuesPerRow,
                        mNumRows, mNumColumns, b.getNumRows(), false, alpha,
                        aIA.get(), aJA.get(), aValues.get(), a.getNumValuesPerRow(),
                        bIA.get(), bJA.get(), bValues.get(), b.getNumValuesPerRow() );
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
    LAMA_LOG_INFO( logger,
                   "this = " << alpha << " * A + " << beta << " * B, with " << "A = " << a << ", B = " << b << ", all are ELL" )

    // TODO: Implement for CUDA
    ContextPtr loc = ContextFactory::getContext( Context::Host );

    LAMA_INTERFACE_FN( matrixAddSizes, loc, ELLUtils, MatrixExpBuild )
    LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )
    LAMA_INTERFACE_FN_T( matrixAdd, loc, ELLUtils, MatrixExp, ValueType )

    LAMA_ASSERT_ERROR( &a != this, "matrixAddMatrix: alias of a with this result matrix" )
    LAMA_ASSERT_ERROR( &b != this, "matrixAddMatrix: alias of b with this result matrix" )

    allocate( a.getNumRows(), a.getNumColumns() );

    LAMA_ASSERT_EQUAL_ERROR( mNumRows, b.getNumRows() )
    LAMA_ASSERT_EQUAL_ERROR( mNumColumns, b.getNumColumns() )

    //mDiagonalProperty = ( mNumRows == mNumColumns );

    {
        ReadAccess<IndexType> aIA( a.getIA(), loc );
        ReadAccess<IndexType> aJA( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );

        ReadAccess<IndexType> bIA( b.getIA(), loc );
        ReadAccess<IndexType> bJA( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );

        WriteOnlyAccess<IndexType> cIA( mIA, loc, mNumRows );

        LAMA_CONTEXT_ACCESS( loc )

        // 1. Step: Compute IA array
        matrixAddSizes( cIA.get(), a.getNumRows(), a.getNumColumns(), false,
                        aIA.get(), aJA.get(), a.getNumValuesPerRow(),
                        bIA.get(), bJA.get(), b.getNumValuesPerRow() );

        // 2. Step: compute length of longest row
        mNumValuesPerRow = maxval( cIA.get(), mNumRows );

        // 3. Step: Allocate IA and Values arrays with new size
        WriteOnlyAccess<IndexType> cJA( mJA, loc, mNumValuesPerRow * mNumRows );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValuesPerRow * mNumRows );

        // 4. Step: Compute cJA and cValues
        matrixAdd( cJA.get(), cValues.get(), cIA.get(), mNumValuesPerRow,
                   mNumRows, mNumColumns, false,
                   alpha, aIA.get(), aJA.get(), aValues.get(), a.getNumValuesPerRow(),
                   beta, bIA.get(), bJA.get(), bValues.get(), b.getNumValuesPerRow() );

    }

    // 5. Step: Computation of C might have produced some zero elements
    compress();

    check( "result of matrix + matrix" ); // just verify for a correct matrix
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>* ELLStorage<ValueType>::create() const
{
    LAMA_LOG_INFO( logger, "create" )

    return new ELLStorage<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>* ELLStorage<ValueType>::copy() const
{
    LAMA_LOG_INFO( logger, "copy" )

    return new ELLStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template class LAMA_DLL_IMPORTEXPORT ELLStorage<float> ;
template class LAMA_DLL_IMPORTEXPORT ELLStorage<double> ;

} // namespace LAMA
