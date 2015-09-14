/**
 * @file JDSStorage.cpp
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
 * @brief Instantitions for template class JDSStorage.
 * @author Thomas Brandes
 * @date 24.06.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/storage/JDSStorage.hpp>

// local library
#include <scai/lama/LAMAInterface.hpp>

#include <scai/lama/LAMAArrayUtils.hpp>

// local scai libraries
#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>

// boost
#include <boost/preprocessor.hpp>

using namespace scai::tasking;
using namespace scai::hmemo;

namespace scai
{

using common::shared_ptr;

namespace lama
{
// Allow for shared_ptr<ValueType> instead of common::shared_ptr<ValueType>


/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSStorage<ValueType>::logger, "MatrixStorage.JDSStorage" )

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const IndexType numRows, const IndexType numColumns )

    : CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>( numRows, numColumns ), mNumDiagonals( 0 ), mNumValues(
          0 )
{
    SCAI_LOG_DEBUG( logger, "JDSStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, IndexType )
    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )

    if( numRows <= 0 )
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

    : CRTPMatrixStorage<JDSStorage<ValueType>,ValueType>( numRows, numColumns ), mNumDiagonals(
          numDiagonals ), mNumValues( numValues ), mDlg( dlg ), mIlg( ilg ), mPerm( perm ), mJa(
          ja ), mValues( values )
{
    check( "JDSStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
    this->resetDiagonalProperty();
    SCAI_LOG_INFO( logger, *this << ": constructed by JDS arrays dlg, ilg, .., values" )
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
    const ContextArray& values )
{
    SCAI_ASSERT_EQUAL_ERROR( numRows, ilg.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows, perm.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numValues, values.size() )
    SCAI_ASSERT_EQUAL_ERROR( numDiagonals, dlg.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumDiagonals = numDiagonals;
    mNumValues = numValues;

    ContextPtr loc = getContextPtr();

    LAMAArrayUtils::assignImpl( mDlg, dlg, loc );
    LAMAArrayUtils::assignImpl( mIlg, ilg, loc );
    LAMAArrayUtils::assignImpl( mPerm, perm, loc );
    LAMAArrayUtils::assignImpl( mJa, ja, loc );

    LAMAArrayUtils::assign( mValues, values, loc ); // supports type conversion

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode

#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "JDSStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif

    this->resetDiagonalProperty();

    SCAI_LOG_INFO( logger, *this << ": set JDS by arrays dlg, ilg, .., values" )
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
    SCAI_LOG_DEBUG( logger, "JDSStorage, matrix is 0 x 0." )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageFormat JDSStorage<ValueType>::getFormat() const
{
    return Format::JDS;
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

    SCAI_LOG_INFO( logger, "setDiagonalImpl with scalar = " << scalar )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( setVal, loc, Utils, Setter, ValueType )

    IndexType numDiagonalValues = std::min( mNumColumns, mNumRows );

    // Note: diagonal is first column in mValues ( stored column-wise )
    // values[i] = scalar

    WriteAccess<ValueType> wValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    setVal( wValues.get(), numDiagonalValues, scalar.getValue<ValueType>() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setDiagonalImpl( const LAMAArray<OtherValueType>& diagonal )
{
    // diagonal property has already been checked

    SCAI_LOG_INFO( logger, "setDiagonalImpl" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( setGather, loc, Utils, Copy, ValueType, OtherValueType )

    IndexType numDiagonal = std::min( mNumColumns, mNumRows );

    ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rJa( mJa, loc );
    WriteOnlyAccess<ValueType> wValues( mValues, loc, numDiagonal );

    SCAI_CONTEXT_ACCESS( loc )

    // diagonal is first column in JDS data
    // values[i] = diagonal[ ja[ i ] ]
    setGather( wValues.get(), rDiagonal.get(), rJa.get(), numDiagonal );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::getRowImpl( LAMAArray<OtherValueType>& row, const IndexType i ) const
{
    SCAI_LOG_INFO( logger, "getRowImpl with i = " << i )

    ContextPtr loc = getContextPtr();

    SCAI_ASSERT_DEBUG( i >= 0 && i < mNumRows, "row index " << i << " out of range" )

    LAMA_INTERFACE_FN_TT( getRow, loc, JDSUtils, Getter, ValueType, OtherValueType )

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );
    ReadAccess<ValueType> values( mValues, loc );
    WriteOnlyAccess<OtherValueType> wRow( row, loc, mNumColumns );

    SCAI_CONTEXT_ACCESS( loc )

    getRow( wRow.get(), i, mNumColumns, mNumRows, perm.get(), ilg.get(), dlg.get(), ja.get(), values.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::getDiagonalImpl( LAMAArray<OtherValueType>& diagonal ) const
{
    SCAI_LOG_INFO( logger, "getDiagonalImpl" )
    //TODO: check diagonal property?
    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( setScatter, loc, Utils, Copy, OtherValueType, ValueType )

    IndexType numDiagonal = std::min( mNumColumns, mNumRows );

    WriteOnlyAccess<OtherValueType> wDiagonal( diagonal, loc, numDiagonal );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<ValueType> rValues( mValues, loc );

    // diagonal is first column in JDS data
    // wDiagonal[ rJa[ i ] ] = rValues[ i ];

    SCAI_CONTEXT_ACCESS( loc )

    setScatter( wDiagonal.get(), rPerm.get(), rValues.get(), numDiagonal );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::scaleImpl( const Scalar scalar )
{
    SCAI_LOG_INFO( logger, "scaleImpl with scalar = " << scalar )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( scale, loc, Utils, Transform, ValueType )

    IndexType size = mValues.size();

    SCAI_ASSERT_EQUAL_DEBUG( size, mNumValues )

    WriteAccess<ValueType> wValues( mValues, loc );
    ValueType value = scalar.getValue<ValueType>();

    SCAI_CONTEXT_ACCESS( loc )

    scale( wValues.get(), value, size );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::scaleImpl( const LAMAArray<OtherValueType>& diagonal )
{
    SCAI_LOG_INFO( logger, "scaleImpl" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_TT( scaleValue, loc, JDSUtils, Scale, ValueType, OtherValueType )

    ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<IndexType> rIlg( mIlg, loc );
    ReadAccess<IndexType> rDlg( mDlg, loc );
    WriteAccess<ValueType> wValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    scaleValue( mNumRows, rPerm.get(), rIlg.get(), rDlg.get(), wValues.get(), rDiagonal.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
bool JDSStorage<ValueType>::checkDiagonalProperty() const
{
    SCAI_LOG_INFO( logger, "checkDiagonalProperty" )

    IndexType n = std::min( mNumRows, mNumColumns );

    bool diagonalProperty = false; // initialization just for safety

    if( n == 0 )
    {
        diagonalProperty = true;
    }
    else if( mNumDiagonals == 0 )
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

        SCAI_CONTEXT_ACCESS( loc )

        diagonalProperty = checkDiagonalProperty( mNumDiagonals, mNumRows, mNumColumns, rPerm.get(), rJa.get(),
                           rDlg.get() );
    }

    return diagonalProperty;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::check( const char* msg ) const
{
    SCAI_LOG_DEBUG( logger, "check at " << getContext() << ", msg = " << msg )

    SCAI_ASSERT_EQUAL_ERROR( mNumRows, mIlg.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumRows, mPerm.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mJa.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumDiagonals, mDlg.size() )

    // check column indexes in JA

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT( validIndexes, loc, Utils, Indexes )

        ReadAccess<IndexType> rJA( mJa, loc );

        SCAI_CONTEXT_ACCESS( loc )

        SCAI_ASSERT_ERROR( validIndexes( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }

    // ToDo: check ILG[0] == mNumDiagonals, be careful about size of ILG

    // check descending values in ILG, DLG

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT_T( isSorted, loc, Utils, Reductions, IndexType )

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        SCAI_CONTEXT_ACCESS( loc )

        bool ascending = false; // check for descending

        SCAI_ASSERT_ERROR( isSorted( rIlg.get(), mNumRows, ascending ),
                           *this << " @ " << msg << ": not descending values in ILG" )

        SCAI_ASSERT_ERROR( isSorted( rDlg.get(), mNumDiagonals, ascending ),
                           *this << " @ " << msg << ": not descending values in DLG" )
    }

    // both, ILG and DLG, must sum up to mNumValues

    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT_T( sum, loc, Utils, Reductions, IndexType )

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        SCAI_CONTEXT_ACCESS( loc )

        SCAI_ASSERT_EQUAL_ERROR( sum( rIlg.get(), mNumRows ), mNumValues )
        SCAI_ASSERT_EQUAL_ERROR( sum( rDlg.get(), mNumDiagonals ), mNumValues )
    }

    // check index values in Perm for out of range

    if( mNumRows > 0 )
    {
        ContextPtr loc = getContextPtr();

        LAMA_INTERFACE_FN_DEFAULT( validIndexes, loc, Utils, Indexes )

        ReadAccess<IndexType> rPerm( mPerm, loc );

        SCAI_CONTEXT_ACCESS( loc )

        SCAI_ASSERT_ERROR( validIndexes( rPerm.get(), mNumRows, mNumRows ),
                           *this << " @ " << msg << ": illegel indexes in Perm" )
    }

    // check perm: no values out of range, but make sure that it is permutation, e.g. [ 0, 0] is illegal

    if( mNumRows > 0 ) // very important as maxval would not work
    {
        ContextPtr loc = getContextPtr();

        // temporary array for inverse permutation, initialize with mNumRows

        LAMAArray<IndexType> invPermArray( mNumRows, mNumRows );

        LAMA_INTERFACE_FN_DEFAULT( setInversePerm, loc, JDSUtils, Sort )

        LAMA_INTERFACE_FN( validIndexes, loc, Utils, Indexes )
        LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )

        ReadAccess<IndexType> rPerm( mPerm, loc );
        WriteAccess<IndexType> wInversePerm( invPermArray, loc );

        SCAI_CONTEXT_ACCESS( loc )

        // set inverse permutation, should overwrite all values 'mNumRows'

        setInversePerm( wInversePerm.get(), rPerm.get(), mNumRows );

        IndexType maxIndex = maxval( wInversePerm.get(), mNumRows );

        SCAI_ASSERT_ERROR( maxIndex < mNumRows, "Perm array does not cover all row indexes, #rows = " << mNumRows );
    }

    // Note: check is not exhaustive, e.g. it does not check for same column index in one row
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_INFO( logger, "set identity values with size = " << size )

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

    SCAI_CONTEXT_ACCESS( loc )

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
    SCAI_LOG_INFO( logger, "setupData" )

    SCAI_ASSERT_EQUAL_ERROR( mIlg.size(), mNumRows )

    LAMA_INTERFACE_FN_T( getValue, loc, Utils, Getter, IndexType )
    LAMA_INTERFACE_FN( ilg2dlg, loc, JDSUtils, Sort )

    ReadAccess<IndexType> ilg( mIlg, loc );
    WriteOnlyAccess<IndexType> dlg( mDlg, loc, mNumDiagonals );

    mNumDiagonals = 0;
    mNumValues = 0;

    SCAI_CONTEXT_ACCESS( loc )

    if( mNumRows )
    {
        mNumDiagonals = getValue( ilg.get(), 0 );
    }

    mNumValues = ilg2dlg( dlg.get(), mNumDiagonals, ilg.get(), mNumRows );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::sortRows( ContextPtr loc )
{
    SCAI_LOG_INFO( logger, *this << "sortRows, number of jagged diagonals = " << mNumDiagonals )

    LAMA_INTERFACE_FN_T( maxval, loc, Utils, Reductions, IndexType )
    LAMA_INTERFACE_FN( sortRows, loc, JDSUtils, Sort )

    // sort the rows according to the array ilg, take sorting over in perm
    WriteAccess<IndexType> ilg( mIlg, loc );
    WriteAccess<IndexType> perm( mPerm, loc );

    SCAI_CONTEXT_ACCESS( loc )

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
    SCAI_REGION( "Storage.JDS->CSR" )

    SCAI_LOG_INFO( logger,
                   "buildCSR<" << common::getScalarType<OtherValueType>() << ">" << " from JDS<" << common::getScalarType<ValueType>() << ">" << " on " << *loc )

    LAMA_INTERFACE_FN_TT( setScatter, loc, Utils, Copy, IndexType, IndexType )
    LAMA_INTERFACE_FN_TT( getCSRValues, loc, JDSUtils, Conversions, ValueType, OtherValueType )
    LAMA_INTERFACE_FN( sizes2offsets, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN( setInversePerm, loc, JDSUtils, Sort )

    ReadAccess<IndexType> rJdsPerm( mPerm, loc );
    ReadAccess<IndexType> rJdsILG( mIlg, loc );
    WriteOnlyAccess<IndexType> wCsrIA( ia, loc, mNumRows + 1 );

    SCAI_CONTEXT_ACCESS( loc )

    // rowValues[ perm[i] ] = ilg[i]
    setScatter( wCsrIA.get(), rJdsPerm.get(), rJdsILG.get(), mNumRows );

    if( ja == NULL || values == NULL )
    {
        wCsrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets( wCsrIA.get(), mNumRows );

    SCAI_ASSERT_EQUAL_DEBUG( numValues, mNumValues )

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
    SCAI_REGION( "Storage.JDS<-CSR" )

    SCAI_LOG_INFO( logger,
                   "setCSRDataImpl<" << common::getScalarType<ValueType>() << "," << common::getScalarType<OtherValueType>() << ">" << ", shape is " << numRows << " x " << numColumns << ", #values for CSR = " << numValues )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN( offsets2sizes, loc, CSRUtils, Offsets )
    LAMA_INTERFACE_FN_T( setOrder, loc, Utils, Setter, IndexType )
    LAMA_INTERFACE_FN_TT( setCSRValues, loc, JDSUtils, Conversions, ValueType, OtherValueType )

    ReadAccess<IndexType> rCsrIA( ia, loc );
    ReadAccess<IndexType> rCsrJA( ja, loc );
    ReadAccess<OtherValueType> rCsrValues( values, loc );

    _MatrixStorage::setDimension( numRows, numColumns );

    mNumValues = numValues;

    // Step 1: fill up the array ilg and perm, detect diagonal property in csr data

    mDiagonalProperty = true; // will be set to false if not valid in one row

    {
        WriteOnlyAccess<IndexType> wIlg( mIlg, loc, mNumRows );
        WriteOnlyAccess<IndexType> wPerm( mPerm, loc, mNumRows );

        SCAI_CONTEXT_ACCESS( loc )

        // ilg willl contain the sizes of each row
        offsets2sizes( wIlg.get(), rCsrIA.get(), mNumRows );

        // set perm to identity permutation
        setOrder( wPerm.get(), mNumRows );
    }

    sortRows( loc ); // sorts ilg and builds perm
    setupData( loc ); // sets dlg, allocates mValues, mJa

    IndexType numDiagonals = mNumDiagonals; // now available

    {
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );

        WriteOnlyAccess<ValueType> wValues( mValues, loc, mNumValues );
        WriteOnlyAccess<IndexType> wJa( mJa, loc, mNumValues );

        SCAI_CONTEXT_ACCESS( loc )

        setCSRValues( wJa.get(), wValues.get(), numRows, rPerm.get(), rIlg.get(), numDiagonals, rDlg.get(),
                      rCsrIA.get(), rCsrJA.get(), rCsrValues.get() );
    }

    this->resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::~JDSStorage()
{
    SCAI_LOG_DEBUG( logger,
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
    SCAI_LOG_INFO( logger, "allocate JDS sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    mNumRows = numRows;
    mNumColumns = numColumns;

    if( mNumRows > 0 )
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

        SCAI_CONTEXT_ACCESS( loc )

        setVal( ilg.get(), mNumRows, 0 );
        setOrder( perm.get(), mNumRows );
    }

    mDiagonalProperty = checkDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "JDSStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << mNumRows << " x " << mNumColumns
           << ", jd = " << mNumDiagonals << ", nnz = " << mNumValues;

    if ( Printable::extended )
    {
        stream << ", context = " << getContext() << ", dlg = " << mDlg << ", ilg = " << mIlg << ", perm = " << mPerm
               << ", ja = " << mJa << ", vales = " << mValues;
    }

    stream << " )";
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ") from " << *this )

    ContextPtr loc = getContextPtr();

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );
    ReadAccess<ValueType> values( mValues, loc );

    LAMA_INTERFACE_FN_TT( getValue, loc, JDSUtils, Getter, ValueType, ValueType )

    SCAI_CONTEXT_ACCESS( loc )

    return getValue( i, j, mNumRows, dlg.get(), ilg.get(), perm.get(), ja.get(), values.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::matrixTimesVector(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.JDS.timesVector" )

    SCAI_LOG_DEBUG( logger,
                    "Computing z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    ContextPtr loc = getContextPtr();

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEMV, loc, JDSUtils, Mult, ValueType )

    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( &result == &y )
    {
        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        // this call will finish the computation, syncToken == NULL

        normalGEMV( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumRows, jdsPerm.get(), jdsILG.get(),
                    mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), NULL );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        // this call will finish the computation, syncToken == NULL

        normalGEMV( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, jdsPerm.get(), jdsILG.get(),
                    mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), NULL );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::vectorTimesMatrix(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.JDS.vectorTimesMatrix" )

    SCAI_LOG_DEBUG( logger,
                    "Computing z = " << alpha << " * x * A + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )

    ContextPtr loc = getContextPtr();

    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEVM, loc, JDSUtils, Mult, ValueType )

    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( &result == &y )
    {
        WriteAccess<ValueType> wResult( result, loc );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        // this call will finish the computation, syncToken == NULL

        normalGEVM( wResult.get(), alpha, rX.get(), beta, wResult.get(), mNumColumns, jdsPerm.get(), jdsILG.get(),
                    mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), NULL );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        ReadAccess<ValueType> rY( y, loc );

        SCAI_CONTEXT_ACCESS( loc )

        // this call will finish the computation, syncToken == NULL

        normalGEVM( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumColumns, jdsPerm.get(), jdsILG.get(),
                    mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), NULL );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
SyncToken* JDSStorage<ValueType>::matrixTimesVectorAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    ContextPtr loc = getContextPtr();

    if ( loc->getType() == context::Host )
    {
        // workaround as common::bind has limited number of arguments and cannot be
        // used later in OpenMP to generate a TaskSyncToken

        void (JDSStorage::*mv)(
            LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>&,
            const ValueType,
            const LAMAArray<ValueType>& ) const

            = &JDSStorage<ValueType>::matrixTimesVector;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        return new TaskSyncToken( bind( mv, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task

    SCAI_REGION( "Storage.JDS.timesVectorAsync" )

    SCAI_LOG_INFO( logger,
                   "Async start z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEMV, loc, JDSUtils, Mult, ValueType )

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > jdsPerm( new ReadAccess<IndexType>( mPerm, loc ) );
    shared_ptr<ReadAccess<IndexType> > jdsDLG( new ReadAccess<IndexType>( mDlg, loc ) );
    shared_ptr<ReadAccess<IndexType> > jdsILG( new ReadAccess<IndexType>( mIlg, loc ) );
    shared_ptr<ReadAccess<IndexType> > jdsJA( new ReadAccess<IndexType>( mJa, loc ) );
    shared_ptr<ReadAccess<ValueType> > jdsValues( new ReadAccess<ValueType>( mValues, loc ) );

    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( &result == &y )
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        syncToken->pushToken( wResult );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        // this call will only start the computation

        normalGEMV( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumRows, jdsPerm->get(), jdsILG->get(),
                    mNumDiagonals, jdsDLG->get(), jdsJA->get(), jdsValues->get(), syncToken.get() );
    }
    else
    {
        shared_ptr<WriteOnlyAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumRows ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        syncToken->pushToken( wResult );
        syncToken->pushToken( rY );

        SCAI_CONTEXT_ACCESS( loc )

        // this call will only start the computation

        normalGEMV( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumRows, jdsPerm->get(), jdsILG->get(),
                    mNumDiagonals, jdsDLG->get(), jdsJA->get(), jdsValues->get(), syncToken.get() );
    }

    syncToken->pushToken( jdsPerm );
    syncToken->pushToken( jdsDLG );
    syncToken->pushToken( jdsILG );
    syncToken->pushToken( jdsJA );
    syncToken->pushToken( jdsValues );
    syncToken->pushToken( rX );

    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
SyncToken* JDSStorage<ValueType>::vectorTimesMatrixAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    ContextPtr loc = getContextPtr();

    if ( loc->getType() == context::Host )
    {
        // workaround as common::bind has limited number of arguments and cannot be
        // used later in OpenMP to generate a TaskSyncToken

        void (JDSStorage::*vm)(
            LAMAArray<ValueType>& result,
            const ValueType alpha,
            const LAMAArray<ValueType>& x,
            const ValueType beta,
            const LAMAArray<ValueType>& y ) const

            = &JDSStorage<ValueType>::vectorTimesMatrix;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        return new TaskSyncToken( bind( vm, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task

    SCAI_REGION( "Storage.JDS.vectorTimesMatrixAsync" )

    SCAI_LOG_INFO( logger,
                   "Async start z = " << alpha << " * x * A + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )

    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )

    LAMA_INTERFACE_FN_T( normalGEVM, loc, JDSUtils, Mult, ValueType )

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.

    shared_ptr<ReadAccess<IndexType> > jdsPerm( new ReadAccess<IndexType>( mPerm, loc ) );
    shared_ptr<ReadAccess<IndexType> > jdsDLG( new ReadAccess<IndexType>( mDlg, loc ) );
    shared_ptr<ReadAccess<IndexType> > jdsILG( new ReadAccess<IndexType>( mIlg, loc ) );
    shared_ptr<ReadAccess<IndexType> > jdsJA( new ReadAccess<IndexType>( mJa, loc ) );
    shared_ptr<ReadAccess<ValueType> > jdsValues( new ReadAccess<ValueType>( mValues, loc ) );

    shared_ptr<ReadAccess<ValueType> > rX( new ReadAccess<ValueType>( x, loc ) );

    // Possible alias of result and y must be handled by coressponding accesses

    if ( &result == &y )
    {
        shared_ptr<WriteAccess<ValueType> > wResult( new WriteAccess<ValueType>( result, loc ) );

        syncToken->pushToken( wResult );

        // we assume that normalGEMV can deal with the alias of result, y

        SCAI_CONTEXT_ACCESS( loc )

        // this call will only start the computation

        normalGEVM( wResult->get(), alpha, rX->get(), beta, wResult->get(), mNumColumns, jdsPerm->get(), jdsILG->get(),
                    mNumDiagonals, jdsDLG->get(), jdsJA->get(), jdsValues->get(), syncToken.get() );
    }
    else
    {
        shared_ptr<WriteOnlyAccess<ValueType> > wResult( new WriteOnlyAccess<ValueType>( result, loc, mNumColumns ) );
        shared_ptr<ReadAccess<ValueType> > rY( new ReadAccess<ValueType>( y, loc ) );

        syncToken->pushToken( wResult );
        syncToken->pushToken( rY );

        SCAI_CONTEXT_ACCESS( loc )

        // this call will only start the computation

        normalGEVM( wResult->get(), alpha, rX->get(), beta, rY->get(), mNumColumns, jdsPerm->get(), jdsILG->get(),
                    mNumDiagonals, jdsDLG->get(), jdsJA->get(), jdsValues->get(), syncToken.get() );
    }

    syncToken->pushToken( jdsPerm );
    syncToken->pushToken( jdsDLG );
    syncToken->pushToken( jdsILG );
    syncToken->pushToken( jdsJA );
    syncToken->pushToken( jdsValues );
    syncToken->pushToken( rX );

    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.JDS.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobi, loc, JDSUtils, Solver, ValueType )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
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

        SCAI_CONTEXT_ACCESS( loc )

        jacobi( wSolution.get(), mNumRows, jdsPerm.get(), jdsIlg.get(), mNumDiagonals, jdsDlg.get(), jdsJA.get(),
                jdsValues.get(), rOldSolution.get(), rRhs.get(), omega, NULL );
    }

}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
SyncToken* JDSStorage<ValueType>::jacobiIterateAsync(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.JDS.jacobiIterateAsync" )

    ContextPtr loc = getContextPtr();

    if( loc->getType() == context::Host )
    {
        // On host we start directly a new task, avoids pushing accesses

        void (JDSStorage::*jb)(
            LAMAArray<ValueType>&,
            const LAMAArray<ValueType>&,
            const LAMAArray<ValueType>&,
            const ValueType omega ) const

            = &JDSStorage<ValueType>::jacobiIterate;

        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;

        return new TaskSyncToken( bind( jb, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    LAMA_INTERFACE_FN_T( jacobi, loc, JDSUtils, Solver, ValueType )

    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, solution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square

    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );

    shared_ptr<WriteAccess<ValueType> > wSolution( new WriteAccess<ValueType>( solution, loc ) );
    syncToken->pushToken( wSolution );
    shared_ptr<ReadAccess<IndexType> > jdsDLG( new ReadAccess<IndexType>( mDlg, loc ) );
    syncToken->pushToken( jdsDLG );
    shared_ptr<ReadAccess<IndexType> > jdsILG( new ReadAccess<IndexType>( mIlg, loc ) );
    syncToken->pushToken( jdsILG );
    shared_ptr<ReadAccess<IndexType> > jdsPerm( new ReadAccess<IndexType>( mPerm, loc ) );
    syncToken->pushToken( jdsPerm );
    shared_ptr<ReadAccess<IndexType> > jdsJA( new ReadAccess<IndexType>( mJa, loc ) );
    syncToken->pushToken( jdsJA );
    shared_ptr<ReadAccess<ValueType> > jdsValues( new ReadAccess<ValueType>( mValues, loc ) );
    syncToken->pushToken( jdsValues );
    shared_ptr<ReadAccess<ValueType> > rOldSolution( new ReadAccess<ValueType>( oldSolution, loc ) );
    syncToken->pushToken( rOldSolution );
    shared_ptr<ReadAccess<ValueType> > rRhs( new ReadAccess<ValueType>( rhs, loc ) );
    syncToken->pushToken( rRhs );

    SCAI_CONTEXT_ACCESS( loc )

    jacobi( wSolution->get(), mNumRows, jdsPerm->get(), jdsILG->get(), mNumDiagonals, jdsDLG->get(), jdsJA->get(),
            jdsValues->get(), rOldSolution->get(), rRhs->get(), omega, syncToken.get() );

    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    SCAI_REGION( "Storage.JDS.jacobiIterateHalo" )

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumRows() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localStorage.getNumColumns() )
    SCAI_ASSERT_DEBUG( localStorage.hasDiagonalProperty(), localStorage << ": has not diagonal property" )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    // need diagonal of local storage in *natural* order
    const LAMAArray<ValueType>* localDiagonal;
    shared_ptr<LAMAArray<ValueType> > tmpLocalDiagonal;
    tmpLocalDiagonal = shared_ptr<LAMAArray<ValueType> >( new LAMAArray<ValueType>() );
    localStorage.getDiagonal( *tmpLocalDiagonal );
    localDiagonal = tmpLocalDiagonal.get();

    jacobiIterateHalo( localSolution, *localDiagonal, oldHaloSolution, omega );
}

/* ------------------------------------------------------------------------------------------------------------------ */
template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const LAMAArray<ValueType>& localDiagonal,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )

    SCAI_REGION( "Storage.JDS.jacobiIterateHalo" )

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( jacobiHalo, loc, JDSUtils, Solver, ValueType )
    LAMA_INTERFACE_FN_T( invert, loc, Utils, Math, ValueType )

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumColumns, oldHaloSolution.size() )

    WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
    ReadAccess<ValueType> diagonal( localDiagonal, loc );
    ReadAccess<IndexType> jdsHaloPerm( mPerm, loc );
    ReadAccess<IndexType> jdsHaloIlg( mIlg, loc );
    ReadAccess<IndexType> jdsHaloDlg( mDlg, loc );
    ReadAccess<IndexType> jdsHaloJA( mJa, loc );
    ReadAccess<ValueType> jdsHaloValues( mValues, loc );
    ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

    SCAI_CONTEXT_ACCESS( loc )

    jacobiHalo( wSolution.get(), mNumRows, diagonal.get(), mNumDiagonals, jdsHaloPerm.get(), jdsHaloIlg.get(),
                jdsHaloDlg.get(), jdsHaloJA.get(), jdsHaloValues.get(), rOldHaloSolution.get(), omega, NULL );

}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType JDSStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    const IndexType n = mNumValues;

    if( n == 0 )
    {
        return 0.0f;
    }

	ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( asum, loc, BLAS, BLAS1, ValueType )

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc )

	return asum( n, data.get(), 1, NULL );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType JDSStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    const IndexType n = mNumValues;

    if( n == 0 )
    {
        return 0.0f;
    }

	ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_T( dot, loc, BLAS, BLAS1, ValueType )

	ReadAccess<ValueType> data( mValues, loc );

	SCAI_CONTEXT_ACCESS( loc )

	return ::sqrt(dot( n, data.get(), 1, data.get(), 1, NULL ));
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    const IndexType n = mNumValues;

    if( n == 0 )
    {
        return 0.0f;
    }

    ContextPtr loc = getContextPtr();

    LAMA_INTERFACE_FN_DEFAULT_T( absMaxVal, loc, Utils, Reductions, ValueType )

    ReadAccess<ValueType> jdsValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

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
    ReadAccess<IndexType> perm( mPerm );
    ReadAccess<IndexType> ilg( mIlg );
    ReadAccess<IndexType> dlg( mDlg );
    ReadAccess<IndexType> ja( mJa );
    ReadAccess<ValueType> values( mValues );

    for( IndexType ii = 0; ii < mNumRows; ii++ )
    {
        cout << "   row " << ii << " is original row " << perm[ii];
        cout << ", #non-zero values = " << ilg[ii] << endl;
        IndexType offset = ii;
        cout << "     column indexes = ";

        for( IndexType d = 0; d < ilg[ii]; d++ )
        {
            cout << " " << ja[offset];
            offset += dlg[d];
        }

        cout << endl;

        offset = ii;
        cout << "     values   = ";

        for( IndexType d = 0; d < ilg[ii]; d++ )
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
JDSStorage<ValueType>* JDSStorage<ValueType>::clone() const
{
    return new JDSStorage();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>* JDSStorage<ValueType>::copy() const
{
    return new JDSStorage( *this );
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

#define LAMA_JDS_STORAGE_INSTANTIATE(z, I, _)                                     \
    template<>                                                                    \
    const char* JDSStorage<ARITHMETIC_HOST_TYPE_##I>::typeName()                  \
    {                                                                             \
        return "JDSStorage<ARITHMETIC_HOST_TYPE_##I>";                            \
    }                                                                             \
                                                                                  \
    template class COMMON_DLL_IMPORTEXPORT JDSStorage<ARITHMETIC_HOST_TYPE_##I> ;

BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_JDS_STORAGE_INSTANTIATE, _ )

#undef LAMA_JDS_STORAGE_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */
