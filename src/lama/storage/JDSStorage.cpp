/**
 * @file JDSStorage.cpp
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
 * @brief Instantitions for template class JDSStorage.
 * @author Thomas Brandes
 * @date 24.06.2011
 * @since 1.0.0
 */

// hpp
#include <lama/storage/JDSStorage.hpp>

// others
#include <lama/ContextFactory.hpp>
#include <lama/ContextAccess.hpp>
#include <lama/LAMAInterface.hpp>

#include <lama/HostReadAccess.hpp>
#include <lama/ReadAccess.hpp>
#include <lama/WriteAccess.hpp>
#include <lama/LAMAArrayUtils.hpp>
#include <lama/task/TaskSyncToken.hpp>

// tracing
#include <lama/tracing.hpp>

namespace lama
{

using boost::shared_ptr;

/* ------------------------------------------------------------------------------------------------------------------ */

LAMA_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSStorage<ValueType>::logger, "MatrixStorage.JDSStorage" )

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const IndexType numRows,
                                   const IndexType numColumns )

    : CRTPMatrixStorage<JDSStorage<ValueType>, ValueType>( numRows, numColumns ),
      mNumDiagonals( 0 ),
      mNumValues( 0 )
{
    LAMA_LOG_DEBUG( logger, "JDSStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )
    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )

    if ( numRows <= 0 )
    {
        return;
    }

    WriteOnlyAccess<IndexType> ilg( mIlg, loc, mNumRows );
    WriteOnlyAccess<IndexType> perm( mPerm, loc, mNumRows );

    setVal( ilg.get(), mNumRows, 0 );
    setOrder( perm.get(), mNumRows );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType numDiagonals,
    const LAMAArray<IndexType>& dlg,
    const LAMAArray<IndexType>& ilg,
    const LAMAArray<IndexType>& perm,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<ValueType>& values )

    : CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>( numRows, numColumns ), 
      mNumDiagonals( numDiagonals ), 
      mNumValues( numValues ), 
      mDlg( dlg ), 
      mIlg( ilg ), 
      mPerm( perm ), 
      mJa( ja ), 
      mValues( values )
{
    check( "JDSStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
    this->resetDiagonalProperty();
    LAMA_LOG_INFO( logger, *this << ": constructed by JDS arrays dlg, ilg, .., values" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setJDSData(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType numDiagonals,
    const LAMAArray<IndexType>& dlg,
    const LAMAArray<IndexType>& ilg,
    const LAMAArray<IndexType>& perm,
    const LAMAArray<IndexType>& ja,
    const _LAMAArray& values )
{
    LAMA_ASSERT_EQUAL_ERROR( numRows, ilg.size() )
    LAMA_ASSERT_EQUAL_ERROR( numRows, perm.size() )
    LAMA_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    LAMA_ASSERT_EQUAL_ERROR( numValues, values.size() )
    LAMA_ASSERT_EQUAL_ERROR( numDiagonals, dlg.size() )

    _MatrixStorage::init( numRows, numColumns );

    mNumDiagonals = numDiagonals;
    mNumValues    = numValues;

    ContextPtr loc = getContextPtr();

    LAMAArrayUtils::assignImpl( mDlg, dlg, loc );
    LAMAArrayUtils::assignImpl( mIlg, ilg, loc );
    LAMAArrayUtils::assignImpl( mPerm, perm, loc );
    LAMAArrayUtils::assignImpl( mJa, ja, loc );

    LAMAArrayUtils::assign( mValues, values, loc );  // supports type conversion

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode

#ifdef LAMA_ASSERT_LEVEL_DEBUG
    check( "JDSStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif

    this->resetDiagonalProperty();

    LAMA_LOG_INFO( logger, *this << ": set JDS by arrays dlg, ilg, .., values" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const JDSStorage<ValueType>& other )

    : CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>( 0, 0 )
{
    assignJDS( other );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const _MatrixStorage& other )

    : CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>( 0, 0 )
{
    assign( other );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>& JDSStorage<ValueType>::operator=( const _MatrixStorage& other )
{
    assign( other );
    return *this;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>& JDSStorage<ValueType>::operator=( const JDSStorage<ValueType>& other )
{
    assignJDS( other );
    return *this;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::assignJDS( const JDSStorage<ValueType>& other )
{
    // this assignment copies directly the arrays

    _MatrixStorage::_assign( other ); // assign member variables of base class

    mNumDiagonals = other.mNumDiagonals;
    mNumValues = other.mNumValues;

    mDlg = other.mDlg;
    mIlg = other.mIlg;
    mPerm = other.mPerm;
    mJa = other.mJa;
    mValues = other.mValues;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;
    mNumDiagonals = 0;
    mNumValues = 0;

    // clear all LAMA arrays used for this storage
    mDlg.clear();
    mIlg.clear();
    mPerm.clear();
    mJa.clear();
    mValues.clear();

    mDiagonalProperty = checkDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage()

    : CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>( 0, 0 ), mNumDiagonals( 0 ), mNumValues( 0 )
{
    LAMA_LOG_DEBUG( logger, "JDSStorage, matrix is 0 x 0." )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageFormat JDSStorage<ValueType>::getFormat() const
{
    return JDS;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType JDSStorage<ValueType>::getNumValues() const
{
    return mNumValues;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType JDSStorage<ValueType>::getNumDiagonals() const
{
    return mNumDiagonals;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setDiagonalImpl( const Scalar scalar )
{
    // diagonal property has already been checked

    LAMA_LOG_INFO( logger, "setDiagonalImpl with scalar = " << scalar )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

    IndexType numDiagonalValues = std::min( mNumColumns, mNumRows );

    // Note: diagonal is first column in mValues ( stored column-wise )
    // values[i] = scalar

    WriteAccess<ValueType> wValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    setVal( wValues.get(), numDiagonalValues, scalar.getValue<ValueType>() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherValueType>& diagonal )
{
    // diagonal property has already been checked

    LAMA_LOG_INFO( logger, "setDiagonalImpl" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( setGather, loc, Utils, Copy, ValueType, OtherValueType )

    IndexType numDiagonal = std::min( mNumColumns, mNumRows );

    ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rJa( mJa, loc );
    WriteOnlyAccess<ValueType> wValues( mValues, loc, numDiagonal );

    LAMA_CONTEXT_ACCESS( loc )

    // diagonal is first column in JDS data
    // values[i] = diagonal[ ja[ i ] ]
    setGather( wValues.get(), rDiagonal.get(), rJa.get(), numDiagonal );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::getRowImpl( LAMAArray<OtherValueType>& row, const IndexType i ) const
{
    LAMA_LOG_INFO( logger, "getRowImpl with i = " << i )

    ContextPtr loc = getContextPtr();

    LAMA_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    LAMA_INTERFACE_FN_TT( getRow, loc, JDSUtils, Getter, ValueType, OtherValueType )

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );
    ReadAccess<ValueType> values( mValues, loc );
    WriteOnlyAccess<OtherValueType> wRow( row, loc, mNumColumns );

    LAMA_CONTEXT_ACCESS( loc )

    getRow( wRow.get(), i, mNumColumns, mNumRows, perm.get(), ilg.get(), dlg.get(), ja.get(), values.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherValueType>& diagonal ) const
{
    LAMA_LOG_INFO( logger, "getDiagonalImpl" )
    //TODO: check diagonal property?
    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( setScatter, loc, Utils, Copy, OtherValueType, ValueType )

    IndexType numDiagonal = std::min( mNumColumns, mNumRows );

    WriteOnlyAccess<OtherValueType> wDiagonal( diagonal, loc, numDiagonal );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<ValueType> rValues( mValues, loc );

    // diagonal is first column in JDS data
    // wDiagonal[ rJa[ i ] ] = rValues[ i ];

    LAMA_CONTEXT_ACCESS( loc )

    setScatter( wDiagonal.get(), rPerm.get(), rValues.get(), numDiagonal );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    LAMA_LOG_INFO( logger, "scaleImpl with scalar = " << scalar )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( scale, loc, Utils, Transform, ValueType, ValueType )

    IndexType size = mValues.size();

    LAMA_ASSERT_EQUAL_DEBUG( size, mNumValues )

    WriteAccess<ValueType> wValues( mValues, loc );
    ValueType value = scalar.getValue<ValueType>();

    LAMA_CONTEXT_ACCESS( loc )

    scale( wValues.get(), size, value );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::scaleImpl( const LAMAArray<OtherValueType>& diagonal )
{
    LAMA_LOG_INFO( logger, "scaleImpl" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( scaleValue, loc, JDSUtils, Scale, ValueType, OtherValueType )

    ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<IndexType> rIlg( mIlg, loc );
    ReadAccess<IndexType> rDlg( mDlg, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    scaleValue( mNumRows, rPerm.get(), rIlg.get(), rDlg.get(), wValues.get(), rDiagonal.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
bool JDSStorage<ValueType>::checkDiagonalProperty() const
{
    LAMA_LOG_INFO( logger, "checkDiagonalProperty" )

    IndexType n = std::min( mNumRows, mNumColumns );

    bool diagonalProperty = false;  // initialization just for safety

    if ( n == 0 ) 
    {
        diagonalProperty = true;
    }
    else if ( mNumDiagonals == 0 )
    {
        // empty storage has no diagonal

        diagonalProperty = false;
    }
    else
    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN( checkDiagonalProperty, loc, JDSUtils, Helper )

        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rJa( mJa, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        LAMA_CONTEXT_ACCESS( loc )

        diagonalProperty = checkDiagonalProperty( mNumDiagonals, mNumRows, mNumColumns, rPerm.get(), rJa.get(), rDlg.get() );
    }

    return diagonalProperty;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::check( const char* msg ) const
{
    LAMA_LOG_DEBUG( logger, "check at " << getContext() << ", msg = " << msg )

    LAMA_ASSERT_EQUAL_ERROR( mNumRows, mIlg.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumRows, mPerm.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mJa.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )
    LAMA_ASSERT_EQUAL_ERROR( mNumDiagonals, mDlg.size() )

    // check column indexes in JA

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT( validIndexes, loc, Utils, Indexes )

        ReadAccess<IndexType> rJA( mJa, loc );

        LAMA_CONTEXT_ACCESS( loc )

        LAMA_ASSERT_ERROR( validIndexes ( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }

    // ToDo: check ILG[0] == mNumDiagonals, be careful about size of ILG

    // check descending values in ILG, DLG

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT_T( isSorted, loc, Utils, Reductions, IndexType )

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        LAMA_CONTEXT_ACCESS( loc )

        bool ascending = false;  // check for descending

        LAMA_ASSERT_ERROR( isSorted ( rIlg.get(), mNumRows, ascending ),
                           *this << " @ " << msg << ": not descending values in ILG" )

        LAMA_ASSERT_ERROR( isSorted ( rDlg.get(), mNumDiagonals, ascending ),
                           *this << " @ " << msg << ": not descending values in DLG" )
    }

    // both, ILG and DLG, must sum up to mNumValues

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT_T( sum, loc, Utils, Reductions, IndexType )

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        LAMA_CONTEXT_ACCESS( loc )

        LAMA_ASSERT_EQUAL_ERROR( sum( rIlg.get(), mNumRows ), mNumValues )
        LAMA_ASSERT_EQUAL_ERROR( sum( rDlg.get(), mNumDiagonals ), mNumValues )
    }

    // check index values in Perm for out of range

    if ( mNumRows > 0 )
    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT( validIndexes, loc, Utils, Indexes )

        ReadAccess<IndexType> rPerm( mPerm, loc );

        LAMA_CONTEXT_ACCESS( loc )

        LAMA_ASSERT_ERROR( validIndexes ( rPerm.get(), mNumRows, mNumRows ),
                           *this << " @ " << msg << ": illegel indexes in Perm" )
    }

    // check perm: no values out of range, but make sure that it is permutation, e.g. [ 0, 0] is illegal

    if ( mNumRows > 0 )   // very important as maxval would not work
    {
        ContextPtr loc = getContextPtr();

        // temporary array for inverse permutation, initialize with mNumRows

        LAMAArray<IndexType> invPermArray( mNumRows, mNumRows );

        LAMA_INTERFACE_FN_DEFAULT( setInversePerm, loc, JDSUtils, Sort )

        LAMA_INTERFACE_FN( validIndexes, loc, Utils, Indexes )
        LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )

        ReadAccess<IndexType> rPerm( mPerm, loc );
        WriteAccess<IndexType> wInversePerm( invPermArray, loc );

        LAMA_CONTEXT_ACCESS( loc )

        // set inverse permutation, should overwrite all values 'mNumRows'

        setInversePerm( wInversePerm.get(), rPerm.get(), mNumRows );

        IndexType maxIndex = maxval( wInversePerm.get(), mNumRows );

        LAMA_ASSERT_ERROR ( maxIndex < mNumRows, "Perm array does not cover all row indexes, #rows = " << mNumRows );
    }

    // Note: check is not exhaustive, e.g. it does not check for same column index in one row
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setIdentity( const IndexType size )
{
    LAMA_LOG_INFO( logger, "set identity values with size = " << size )

    ContextPtr loc = getContextPtr();

    mNumRows = size;
    mNumColumns = size;
    mNumDiagonals = 1; // identity has exactly one diagonal
    mNumValues = mNumRows;

    WriteOnlyAccess<IndexType> wDlg( mDlg, loc, mNumDiagonals );
    WriteOnlyAccess<IndexType> wIlg( mIlg, loc, mNumRows );
    WriteOnlyAccess<IndexType> wPerm( mPerm, loc, mNumRows );
    WriteOnlyAccess<IndexType> wJa( mJa, loc, mNumValues );
    WriteOnlyAccess<ValueType> wValues( mValues, loc, mNumValues );

    LAMA_CONTEXT_ACCESS( loc )

    {
        ValueType one = static_cast<ValueType>( 1.0 );
        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

        setVal( wValues.get(), mNumRows, one );
    }

    IndexType one = 1;
    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )
    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )

    setVal( wDlg.get(), 1, mNumRows );
    setVal( wIlg.get(), mNumRows, one );
    setOrder( wPerm.get(), mNumRows );
    setOrder( wJa.get(), mNumRows );

    mDiagonalProperty = true;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setupData( ContextPtr loc )
{
    LAMA_LOG_INFO( logger, "setupData" )

    LAMA_ASSERT_EQUAL_ERROR( mIlg.size(), mNumRows )

    LAMA_INTERFACE_FN_T( getValue, loc, Utils, Getter, IndexType )
    LAMA_INTERFACE_FN( ilg2dlg, loc, JDSUtils, Sort )

    ReadAccess<IndexType> ilg( mIlg, loc );
    WriteOnlyAccess<IndexType> dlg( mDlg, loc, mNumDiagonals );

    mNumDiagonals = 0;
    mNumValues = 0;

    LAMA_CONTEXT_ACCESS( loc )

    if ( mNumRows )
    {
        mNumDiagonals = getValue( ilg.get(), 0 );
    }

    mNumValues = ilg2dlg( dlg.get(), mNumDiagonals, ilg.get(), mNumRows );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::sortRows( ContextPtr loc )
{
    LAMA_LOG_INFO( logger, *this << "sortRows, number of jagged diagonals = " << mNumDiagonals )

    LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )
    LAMA_INTERFACE_FN( sortRows, loc, JDSUtils, Sort )

    // sort the rows according to the array ilg, take sorting over in perm
    WriteAccess<IndexType> ilg( mIlg, loc );
    WriteAccess<IndexType> perm( mPerm, loc );

    LAMA_CONTEXT_ACCESS( loc )

    mNumDiagonals = maxval( ilg.get(), mNumRows );

    sortRows( ilg.get(), perm.get(), mNumRows );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::buildCSR(
    LAMAArray<IndexType>& ia,
    LAMAArray<IndexType>* ja,
    LAMAArray<OtherValueType>* values,
    const ContextPtr loc ) const
{
    LAMA_REGION( "Storage.JDS->CSR" )

    LAMA_LOG_INFO( logger,
                   "buildCSR<" << Scalar::getType<OtherValueType>() << ">" << " from JDS<" << Scalar::getType<ValueType>() << ">" << " on " << *loc )

    LAMA_INTERFACE_FN_TT( setScatter, loc, Utils, Copy, IndexType, IndexType )
    LAMA_INTERFACE_FN_TT( getCSRValues, loc, JDSUtils, Conversions, ValueType, OtherValueType )
    LAMA_INTERFACE_FN( sizes2offsets, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN( setInversePerm, loc, JDSUtils, Sort )

    ReadAccess<IndexType> rJdsPerm( mPerm, loc );
    ReadAccess<IndexType> rJdsILG( mIlg, loc );
    WriteOnlyAccess<IndexType> wCsrIA( ia, loc, mNumRows + 1 );

    LAMA_CONTEXT_ACCESS( loc )

    // rowValues[ perm[i] ] = ilg[i]
    setScatter( wCsrIA.get(), rJdsPerm.get(), rJdsILG.get(), mNumRows );

    if ( ja == NULL || values == NULL )
    {
        wCsrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets( wCsrIA.get(), mNumRows );

    LAMA_ASSERT_EQUAL_DEBUG( numValues, mNumValues )

    // temporary array for inverse permutation
    LAMAArray<IndexType> invPermArray; // allows to find a CSR row in JDS rows

    WriteOnlyAccess<IndexType> wJdsInversePerm( invPermArray, loc, mNumRows );

    // compute the inverse permutation so that we find original row in JDS data
    setInversePerm( wJdsInversePerm.get(), rJdsPerm.get(), mNumRows );

    WriteOnlyAccess<IndexType> wCsrJA( *ja, loc, mNumValues );
    WriteOnlyAccess<OtherValueType> wCsrValues( *values, loc, mNumValues );

    ReadAccess<IndexType> rJdsDLG( mDlg, loc );
    ReadAccess<IndexType> rJdsJA( mJa, loc );
    ReadAccess<ValueType> rJdsValues( mValues, loc );

    // now we can convert JDS to CSR via interface
    getCSRValues( wCsrJA.get(), wCsrValues.get(), wCsrIA.get(), mNumRows, wJdsInversePerm.get(), rJdsILG.get(),
                  rJdsDLG.get(), rJdsJA.get(), rJdsValues.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const LAMAArray<IndexType>& ia,
    const LAMAArray<IndexType>& ja,
    const LAMAArray<OtherValueType>& values,
    const ContextPtr )
{
    LAMA_REGION( "Storage.JDS<-CSR" )

    LAMA_LOG_INFO( logger,
                   "setCSRDataImpl<" << typeid(ValueType).name() << "," << typeid(OtherValueType).name() << ">" << ", shape is " << numRows << " x " << numColumns << ", #values for CSR = " << numValues )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN( offsets2sizes, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )
    LAMA_INTERFACE_FN_TT( setCSRValues, loc, JDSUtils, Conversions, ValueType, OtherValueType )

    ReadAccess<IndexType> rCsrIA( ia, loc );
    ReadAccess<IndexType> rCsrJA( ja, loc );
    ReadAccess<OtherValueType> rCsrValues( values, loc );

    _MatrixStorage::init( numRows, numColumns );

    mNumValues = numValues;

    // Step 1: fill up the array ilg and perm, detect diagonal property in csr data

    mDiagonalProperty = true; // will be set to false if not valid in one row

    {
        WriteOnlyAccess<IndexType> wIlg( mIlg, loc, mNumRows );
        WriteOnlyAccess<IndexType> wPerm( mPerm, loc, mNumRows );

        LAMA_CONTEXT_ACCESS( loc )

        // ilg willl contain the sizes of each row
        offsets2sizes( wIlg.get(), rCsrIA.get(), mNumRows );

        // set perm to identity permutation
        setOrder( wPerm.get(), mNumRows );
    }

    sortRows( loc ); // sorts ilg and builds perm
    setupData( loc ); // sets dlg, allocates mValues, mJa

    IndexType numDiagonals = mNumDiagonals;  // now available

    {
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        WriteOnlyAccess<ValueType> wValues( mValues, loc, mNumValues );
        WriteOnlyAccess<IndexType> wJa( mJa, loc, mNumValues );

        LAMA_CONTEXT_ACCESS( loc )

        setCSRValues( wJa.get(), wValues.get(), numRows, rPerm.get(), rIlg.get(), 
                      numDiagonals, rDlg.get(), rCsrIA.get(), rCsrJA.get(), rCsrValues.get() );
    }

    this->resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::~JDSStorage()
{
    LAMA_LOG_DEBUG( logger,
                    "~JDSStorage for matrix " << mNumRows << " x " << mNumColumns << ", # diags = " << mNumDiagonals )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::purge()
{
    mNumColumns = 0;
    mNumRows = 0;
    mNumDiagonals = 0;
    mNumValues = 0;
    mDiagonalProperty = checkDiagonalProperty();

    mDlg.purge();
    mIlg.purge();
    mPerm.purge();
    mJa.purge();
    mValues.purge();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    LAMA_LOG_INFO( logger, "allocate JDS sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    mNumRows    = numRows;
    mNumColumns = numColumns;

    if ( mNumRows >  0 )
    {
        // the arrays mIlg and mPerm need initialization

        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )
        LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )

        mNumRows = numRows;
        mNumColumns = numColumns;

        // we allocate at least ilg and perm with the correct value for a zero matrix

        WriteOnlyAccess<IndexType> ilg( mIlg, loc, mNumRows );
        WriteOnlyAccess<IndexType> perm( mPerm, loc, mNumRows );

        LAMA_CONTEXT_ACCESS( loc )

        setVal( ilg.get(), mNumRows, 0 );
        setOrder( perm.get(), mNumRows );
    }

    mDiagonalProperty = checkDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "JDS( rows = " << mNumRows << ", cols = " << mNumColumns << ", jd = " << mNumDiagonals << ", values = "
           << mNumValues;

    if ( Printable::extended )
    {
        stream << ", context = " << getContext() << ", dlg = " << mDlg << ", ilg = " << mIlg << ", perm = " << mPerm
               << ", ja = " << mJa << ", vales = " << mValues;
    }

    stream << std::endl;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    LAMA_LOG_TRACE( logger, "get value (" << i << ", " << j << ") from " << *this )

    ContextPtr loc = getContextPtr();

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );
    ReadAccess<ValueType> values( mValues, loc );

    LAMA_INTERFACE_FN_TT( getValue, loc, JDSUtils, Getter, ValueType, ValueType )

    LAMA_CONTEXT_ACCESS( loc )

    return getValue( i, j, mNumRows, dlg.get(), ilg.get(), perm.get(), ja.get(), values.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::matrixTimesVector(

    LAMAArrayView<ValueType> result,
    const ValueType alpha,
    const LAMAArrayConstView<ValueType> x,
    const ValueType beta,
    const LAMAArrayConstView<ValueType> y ) const

{
    // TODO: check CUDA implementation

    LAMA_REGION( "Storage.JDS.timesVector" )

    LAMA_LOG_DEBUG( logger,
                    "Computing z = alpha * A * x + beta * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    LAMA_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    LAMA_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    ContextPtr loc = getContextPtr();

    LAMA_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEMV, loc, JDSUtils, Mult, ValueType )

    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( result == y )
    {
        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        LAMA_CONTEXT_ACCESS( loc )

        normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, jdsPerm.get(), jdsILG.get(),
                    mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), NULL );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        LAMA_CONTEXT_ACCESS( loc )

        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, jdsPerm.get(), jdsILG.get(),
                    mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), NULL );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterate(
    LAMAArrayView<ValueType> solution,
    const LAMAArrayConstView<ValueType> oldSolution,
    const LAMAArrayConstView<ValueType> rhs,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.JDS.jacobiIterate" )

    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobi, loc, JDSUtils, Solver, ValueType )

    LAMA_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( solution == oldSolution )
    {
        LAMA_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square

    {

        WriteAccess<ValueType> wSolution( solution, loc );
        ReadAccess<IndexType> jdsDlg( mDlg, loc );
        ReadAccess<IndexType> jdsIlg( mIlg, loc );
        ReadAccess<IndexType> jdsPerm( mPerm, loc );
        ReadAccess<IndexType> jdsJA( mJa, loc );
        ReadAccess<ValueType> jdsValues( mValues, loc );
        ReadAccess<ValueType> rOldSolution( oldSolution, loc );
        ReadAccess<ValueType> rRhs( rhs, loc );

        LAMA_CONTEXT_ACCESS( loc )

        jacobi( wSolution.get(), mNumRows, jdsPerm.get(), jdsIlg.get(), mNumDiagonals, jdsDlg.get(), jdsJA.get(),
                jdsValues.get(), rOldSolution.get(), rRhs.get(), omega, NULL );
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
SyncToken* JDSStorage<ValueType>::jacobiIterateAsync(
    LAMAArrayView<ValueType> solution,
    const LAMAArrayConstView<ValueType> oldSolution,
    const LAMAArrayConstView<ValueType> rhs,
    const ValueType omega ) const
{
    LAMA_REGION( "Storage.JDS.jacobiIterateAsync" )

    ContextPtr loc = getContextPtr();

    if ( loc->getType() == Context::Host )
    {
        // On host we start directly a new task, avoids pushing accesses

        void ( JDSStorage::*jb )(
            LAMAArrayView<ValueType>,
            const LAMAArrayConstView<ValueType>,
            const LAMAArrayConstView<ValueType>,
            const ValueType omega ) const

        = &JDSStorage<ValueType>::jacobiIterate;

        return new TaskSyncToken( boost::bind( jb, this, solution, oldSolution, rhs, omega ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task

    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    LAMA_INTERFACE_FN_T( jacobi, loc, JDSUtils, Solver, ValueType )

    LAMA_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( solution == oldSolution )
    {
        LAMA_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square

    std::auto_ptr<SyncToken> syncToken( loc->getSyncToken() );

    shared_ptr<WriteAccess<ValueType> > wSolution( new WriteAccess<ValueType>( solution, loc ) );
    syncToken->pushAccess( wSolution );
    shared_ptr<ReadAccess<IndexType> > jdsDLG( new ReadAccess<IndexType>( mDlg, loc ) );
    syncToken->pushAccess( jdsDLG );
    shared_ptr<ReadAccess<IndexType> > jdsILG( new ReadAccess<IndexType>( mIlg, loc ) );
    syncToken->pushAccess( jdsILG );
    shared_ptr<ReadAccess<IndexType> > jdsPerm( new ReadAccess<IndexType>( mPerm, loc ) );
    syncToken->pushAccess( jdsPerm );
    shared_ptr<ReadAccess<IndexType> > jdsJA( new ReadAccess<IndexType>( mJa, loc ) );
    syncToken->pushAccess( jdsJA );
    shared_ptr<ReadAccess<ValueType> > jdsValues( new ReadAccess<ValueType>( mValues, loc ) );
    syncToken->pushAccess( jdsValues );
    shared_ptr<ReadAccess<ValueType> > rOldSolution( new ReadAccess<ValueType>( oldSolution, loc ) );
    syncToken->pushAccess( rOldSolution );
    shared_ptr<ReadAccess<ValueType> > rRhs( new ReadAccess<ValueType>( rhs, loc ) );
    syncToken->pushAccess( rRhs );

    LAMA_CONTEXT_ACCESS( loc )

    jacobi( wSolution->get(), mNumRows, jdsPerm->get(), jdsILG->get(), mNumDiagonals, jdsDLG->get(), jdsJA->get(),
            jdsValues->get(), rOldSolution->get(), rRhs->get(), omega, syncToken.get() );

    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterateHalo(
    LAMAArrayView<ValueType> localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const LAMAArrayConstView<ValueType> oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    LAMA_REGION( "Storage.JDS.jacobiIterateHalo" )

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumRows() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumColumns() )
    LAMA_ASSERT_DEBUG( localStorage.hasDiagonalProperty(), localStorage << ": has not diagonal property" )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    // need diagonal of local storage in *natural* order
    const LAMAArray<ValueType>* localDiagonal;
    boost::shared_ptr<LAMAArray<ValueType> > tmpLocalDiagonal;
    tmpLocalDiagonal = boost::shared_ptr<LAMAArray<ValueType> >( new LAMAArray<ValueType>() );
    localStorage.getDiagonal( *tmpLocalDiagonal );
    localDiagonal = tmpLocalDiagonal.get();

    jacobiIterateHalo( localSolution, localDiagonal, oldHaloSolution, omega );

}

/* ------------------------------------------------------------------------------------------------------------------ */
template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterateHalo(
    LAMAArrayView<ValueType> localSolution,
    const LAMAArray<ValueType>* localDiagonal,
    const LAMAArrayConstView<ValueType> oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    LAMA_REGION( "Storage.JDS.jacobiIterateHalo" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobiHalo, loc, JDSUtils, Solver, ValueType )
    LAMA_INTERFACE_FN_T( invert, loc, Utils, Math, ValueType )

    LAMA_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    LAMA_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
    ReadAccess<ValueType> diagonal( *localDiagonal, loc );
    ReadAccess<IndexType> jdsHaloPerm( mPerm, loc );
    ReadAccess<IndexType> jdsHaloIlg( mIlg, loc );
    ReadAccess<IndexType> jdsHaloDlg( mDlg, loc );
    ReadAccess<IndexType> jdsHaloJA( mJa, loc );
    ReadAccess<ValueType> jdsHaloValues( mValues, loc );
    ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

    LAMA_CONTEXT_ACCESS( loc )

    jacobiHalo( wSolution.get(), mNumRows, diagonal.get(), mNumDiagonals, jdsHaloPerm.get(), jdsHaloIlg.get(),
                jdsHaloDlg.get(), jdsHaloJA.get(), jdsHaloValues.get(), rOldHaloSolution.get(), omega, NULL );

}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::maxNorm() const
{
    LAMA_LOG_INFO( logger, *this << ": maxNorm()" )

    const IndexType n = mNumValues;

    if ( n == 0 )
    {
        return 0.0f;
    }

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, Utils, Reductions, ValueType )

    ReadAccess<ValueType> jdsValues( mValues, loc );

    LAMA_CONTEXT_ACCESS( loc )

    ValueType maxval = absMaxVal( jdsValues.get(), n );

    return maxval;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::print() const
{
    using std::cout;
    using std::endl;

    cout << "JDSStorage of matrix " << mNumRows << " x " << mNumColumns;
    cout << ", #non-zero values = " << mNumValues << endl;
    HostReadAccess<IndexType> perm( mPerm );
    HostReadAccess<IndexType> ilg( mIlg );
    HostReadAccess<IndexType> dlg( mDlg );
    HostReadAccess<IndexType> ja( mJa );
    HostReadAccess<ValueType> values( mValues );

    for ( IndexType ii = 0; ii < mNumRows; ii++ )
    {
        cout << "   row " << ii << " is original row " << perm[ii];
        cout << ", #non-zero values = " << ilg[ii] << endl;
        IndexType offset = ii;
        cout << "     column indexes = ";

        for ( IndexType d = 0; d < ilg[ii]; d++ )
        {
            cout << " " << ja[offset];
            offset += dlg[d];
        }
        cout << endl;

        offset = ii;
        cout << "     values   = ";

        for ( IndexType d = 0; d < ilg[ii]; d++ )
        {
            cout << " " << values[offset];
            offset += dlg[d];
        }

        cout << endl;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mDlg.prefetch( location );
    mIlg.prefetch( location );
    mPerm.prefetch( location );
    mJa.prefetch( location );
    mValues.prefetch( location );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const LAMAArray<IndexType>& JDSStorage<ValueType>::getJA() const
{
    return mJa;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const LAMAArray<IndexType>& JDSStorage<ValueType>::getPerm() const
{
    return mPerm;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const LAMAArray<IndexType>& JDSStorage<ValueType>::getDlg() const
{
    return mDlg;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const LAMAArray<IndexType>& JDSStorage<ValueType>::getIlg() const
{
    return mIlg;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const LAMAArray<ValueType>& JDSStorage<ValueType>::getValues() const
{
    return mValues;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::wait() const
{
    mDlg.wait();
    mIlg.wait();
    mPerm.wait();
    mJa.wait();
    mValues.wait();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::swap( JDSStorage<ValueType>& other )
{
    MatrixStorage<ValueType>::swap( other ); // swap member variable of base class

    std::swap( mNumDiagonals, other.mNumDiagonals );
    std::swap( mNumValues, other.mNumValues );
    mDlg.swap( other.mDlg );
    mIlg.swap( other.mIlg );
    mPerm.swap( other.mPerm );
    mJa.swap( other.mJa );
    mValues.swap( other.mValues );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
size_t JDSStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType);
    memoryUsage += sizeof(IndexType) * mDlg.size();
    memoryUsage += sizeof(IndexType) * mIlg.size();
    memoryUsage += sizeof(IndexType) * mPerm.size();
    memoryUsage += sizeof(IndexType) * mJa.size();
    memoryUsage += sizeof(ValueType) * mValues.size();
    return memoryUsage;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>* JDSStorage<ValueType>::create() const
{
    return new JDSStorage();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>* JDSStorage<ValueType>::copy() const
{
    return new JDSStorage( *this );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<>
const char* JDSStorage<float>::typeName()
{
    return "JDSStorage<float>";
}

/* ------------------------------------------------------------------------------------------------------------------ */

template class LAMA_DLL_IMPORTEXPORT JDSStorage<float> ;

template<>
const char* JDSStorage<double>::typeName()
{
    return "JDSStorage<double>";
}

template class LAMA_DLL_IMPORTEXPORT JDSStorage<double> ;

}
