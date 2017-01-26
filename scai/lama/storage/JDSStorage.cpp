/**
 * @file JDSStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Instantitions for template class JDSStorage.
 * @author Thomas Brandes
 * @date 24.06.2011
 */

// hpp
#include <scai/lama/storage/JDSStorage.hpp>

// local scai libraries
#include <scai/sparsekernel/JDSKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::shared_ptr;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;

using sparsekernel::CSRKernelTrait;
using sparsekernel::JDSKernelTrait;

namespace lama
{
// Allow for shared_ptr<ValueType> instead of common::shared_ptr<ValueType>


/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSStorage<ValueType>::logger, "MatrixStorage.JDSStorage" )

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const IndexType numRows, const IndexType numColumns )

    : CRTPMatrixStorage<JDSStorage<ValueType>, ValueType>( numRows, numColumns ),
      mNumDiagonals( 0 ),
      mNumValues( 0 )
{
    SCAI_LOG_DEBUG( logger, "JDSStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )

    if ( numRows <= 0 )
    {
        return;
    }

    ContextPtr prefLoc = this->getContextPtr();
    mIlg.clear();
    mIlg.resize( mNumRows );
    HArrayUtils::setScalar( mIlg, IndexType( 0 ), utilskernel::binary::COPY, prefLoc );
    HArrayUtils::setOrder( mPerm, mNumRows, prefLoc );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const IndexType numDiagonals,
    const HArray<IndexType>& dlg,
    const HArray<IndexType>& ilg,
    const HArray<IndexType>& perm,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values )

    : CRTPMatrixStorage<JDSStorage<ValueType>, ValueType>( numRows, numColumns ), mNumDiagonals(
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
    const HArray<IndexType>& dlg,
    const HArray<IndexType>& ilg,
    const HArray<IndexType>& perm,
    const HArray<IndexType>& ja,
    const _HArray& values )
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
    HArrayUtils::setArray( mDlg, dlg, utilskernel::binary::COPY, loc );
    HArrayUtils::setArray( mIlg, ilg, utilskernel::binary::COPY, loc );
    HArrayUtils::setArray( mPerm, perm, utilskernel::binary::COPY, loc );
    HArrayUtils::setArray( mJa, ja, utilskernel::binary::COPY, loc );
    HArrayUtils::assign( mValues, values, loc ); // supports type conversion
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

    : CRTPMatrixStorage<JDSStorage<ValueType>, ValueType>( 0, 0 )
{
    assignJDS( other );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const _MatrixStorage& other )

    : CRTPMatrixStorage<JDSStorage<ValueType>, ValueType>( 0, 0 )
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
    mDiagonalProperty = other.mDiagonalProperty;
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

    : CRTPMatrixStorage<JDSStorage<ValueType>, ValueType>( 0, 0 ), mNumDiagonals( 0 ), mNumValues( 0 )
{
    SCAI_LOG_DEBUG( logger, "JDSStorage, matrix is 0 x 0." )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
Format::MatrixStorageFormat JDSStorage<ValueType>::getFormat() const
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
void JDSStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    // diagonal property has already been checked
    SCAI_LOG_INFO( logger, "setDiagonalImpl with value = " << value )
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    ContextPtr loc = this->getContextPtr();
    setVal.getSupportedContext( loc );
    IndexType numDiagonalValues = common::Math::min( mNumColumns, mNumRows );
    // Note: diagonal is first column in mValues ( stored column-wise )
    // values[i] = scalar
    WriteAccess<ValueType> wValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    setVal[loc]( wValues.get(), numDiagonalValues, value, utilskernel::binary::COPY );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setDiagonalImpl( const HArray<OtherValueType>& diagonal )
{
    // diagonal property has already been checked
    SCAI_LOG_INFO( logger, "setDiagonalImpl" )
    static LAMAKernel<UtilKernelTrait::setGather<ValueType, OtherValueType> > setGather;
    ContextPtr loc = this->getContextPtr();
    setGather.getSupportedContext( loc );
    IndexType numDiagonal = common::Math::min( mNumColumns, mNumRows );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rJa( mJa, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    // diagonal is first column in JDS data
    // values[i] = diagonal[ ja[ i ] ]
    setGather[loc]( wValues.get(), rDiagonal.get(), rJa.get(), utilskernel::binary::COPY, numDiagonal );
    // Still problem to use HArrayUtils::gather, as only part of the array is used
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::getRowImpl( HArray<OtherValueType>& row, const IndexType i ) const
{
    SCAI_REGION( "Storage.JDS.getRow" )

    SCAI_LOG_INFO( logger, "getRowImpl with i = " << i )
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    static LAMAKernel<JDSKernelTrait::getRow<ValueType, OtherValueType> > getRow;
    ContextPtr loc = this->getContextPtr();
    getRow.getSupportedContext( loc );
    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );
    ReadAccess<ValueType> values( mValues, loc );
    WriteOnlyAccess<OtherValueType> wRow( row, loc, mNumColumns );
    SCAI_CONTEXT_ACCESS( loc )
    getRow[loc]( wRow.get(), i, mNumColumns, mNumRows, perm.get(), ilg.get(), dlg.get(), ja.get(), values.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::_HArray& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.JDS.getSparseRow" )

    SCAI_LOG_INFO( logger, "getSparseRow( " << i << " of " << mNumRows << " )" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )

    static LAMAKernel<JDSKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = this->getContextPtr();

    getValuePosRow.getSupportedContext( loc );

    HArray<IndexType> pos;  // positions in the array mJa, mValues

    {
        SCAI_CONTEXT_ACCESS( loc )

        // start with maximal possible size, is number of columns, resize later

        WriteOnlyAccess<IndexType> wPos( pos, loc, mNumColumns );  

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );

        IndexType cnt = getValuePosRow[loc]( wPos.get(), i, mNumRows,
                                             rIlg.get(), rDlg.get(), rPerm.get() );

        wPos.resize( cnt );
    }

    // with entries in pos we can gather the column indexes and the values from jdsJA, jdsValues
 
    HArrayUtils::gatherImpl( jA, mJa, pos, utilskernel::binary::COPY, loc );
    HArrayUtils::gather( values, mValues, pos, utilskernel::binary::COPY, loc );

    SCAI_LOG_DEBUG( logger, "getSparseRow( " << i << " ) : jA = " << jA << ", values = " << values )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void JDSStorage<ValueType>::getColumnImpl( HArray<OtherType>& column, const IndexType j ) const
{
    SCAI_LOG_INFO( logger, "getColumn( " << j << " ) of : " << *this )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_REGION( "Storage.JDS.getCol" )

    static LAMAKernel<JDSKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, mNumRows );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumRows );

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rJa( mJa, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j, mNumRows,
                                             rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    SCAI_LOG_INFO( logger, "getColumn( " << j << " ) with " << valuePos.size() << " non-zero entries" )

    column.init( ValueType( 0 ), mNumRows );

    // column[ row ] = mValues[ pos ];

    HArrayUtils::gatherImpl( colValues, mValues, valuePos, utilskernel::binary::COPY, loc );
    HArrayUtils::scatterImpl( column, rowIndexes, true, colValues, utilskernel::binary::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void JDSStorage<ValueType>::setRowImpl( const HArray<OtherType>& row, const IndexType i,
                                        const utilskernel::binary::BinaryOp op )
{
    SCAI_REGION( "Storage.JDS.setRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), mNumColumns, "row array to small for set" )

    SCAI_LOG_INFO( logger, "setRowImpl( i = " << i << " )" )

    static LAMAKernel<JDSKernelTrait::setRow<ValueType, OtherType> > setRow;

    ContextPtr loc = this->getContextPtr();
    setRow.getSupportedContext( loc );

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );
    WriteAccess<ValueType> values( mValues, loc );
    ReadAccess<OtherType> rRow( row, loc );
    SCAI_CONTEXT_ACCESS( loc )
    setRow[loc]( values.get(), i, mNumColumns, mNumRows, perm.get(), ilg.get(), dlg.get(), ja.get(), rRow.get(), op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void JDSStorage<ValueType>::setColumnImpl( const HArray<OtherType>& column, const IndexType j,
        const utilskernel::binary::BinaryOp op )
{
    SCAI_LOG_INFO( logger, "setColumn( " << j << " ) of : " << *this << " with column " << column )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), mNumRows, "column array to small for set" )

    SCAI_REGION( "Storage.JDS.setCol" )

    static LAMAKernel<JDSKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, mNumRows );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumRows );

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rJa( mJa, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j, mNumRows,
                                             rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    SCAI_LOG_INFO( logger, "setColumn( " << j << " ) updates " << rowIndexes.size() << " entries" )

    //  mValues[ pos ] op= column[row]

    HArrayUtils::gatherImpl( colValues, column, rowIndexes, utilskernel::binary::COPY, loc );
    HArrayUtils::scatterImpl( mValues, valuePos, true, colValues, op, loc );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::getDiagonalImpl( HArray<OtherValueType>& diagonal ) const
{
    SCAI_LOG_INFO( logger, "getDiagonalImpl" )
    //TODO: check diagonal property?
    static LAMAKernel<UtilKernelTrait::setScatter<OtherValueType, ValueType> > setScatter;
    ContextPtr loc = this->getContextPtr();
    setScatter.getSupportedContext( loc );
    IndexType numDiagonal = common::Math::min( mNumColumns, mNumRows );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<OtherValueType> wDiagonal( diagonal, loc, numDiagonal );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<ValueType> rValues( mValues, loc );
    // diagonal is first column in JDS data
    // wDiagonal[ rJa[ i ] ] = rValues[ i ];
    setScatter[loc]( wDiagonal.get(), rPerm.get(), true, rValues.get(), utilskernel::binary::COPY, numDiagonal );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::scaleImpl( const ValueType value )
{
    SCAI_LOG_INFO( logger, "scaleImpl with value = " << value )
    HArrayUtils::binaryOpScalar2( mValues, mValues, value, utilskernel::binary::MULT, this->getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */


template<typename ValueType>
void JDSStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, utilskernel::unary::CONJ, this->getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::scaleImpl( const HArray<OtherValueType>& diagonal )
{
    SCAI_LOG_INFO( logger, "scaleImpl" )
    static LAMAKernel<JDSKernelTrait::scaleRows<ValueType, OtherValueType> > scaleRows;
    ContextPtr loc = this->getContextPtr();
    scaleRows.getSupportedContext( loc );
    ReadAccess<OtherValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<IndexType> rIlg( mIlg, loc );
    ReadAccess<IndexType> rDlg( mDlg, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    scaleRows[loc]( wValues.get(), mNumRows, rPerm.get(), rIlg.get(), rDlg.get(), rDiagonal.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
bool JDSStorage<ValueType>::checkDiagonalProperty() const
{
    SCAI_LOG_INFO( logger, "checkDiagonalProperty" )
    IndexType n = common::Math::min( mNumRows, mNumColumns );
    bool diagonalProperty = false; // initialization just for safety

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
        static LAMAKernel<JDSKernelTrait::checkDiagonalProperty> checkDiagonalProperty;
        ContextPtr loc = this->getContextPtr();
        checkDiagonalProperty.getSupportedContext( loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rJa( mJa, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        SCAI_CONTEXT_ACCESS( loc )
        diagonalProperty = checkDiagonalProperty[loc]( mNumDiagonals, mNumRows, mNumColumns, rPerm.get(), rJa.get(),
                           rDlg.get() );
    }

    return diagonalProperty;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::check( const char* msg ) const
{
    SCAI_LOG_DEBUG( logger, "check at " << *getContextPtr() << ", msg = " << msg )
    SCAI_ASSERT_EQUAL_ERROR( mNumRows, mIlg.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumRows, mPerm.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mJa.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumValues, mValues.size() )
    SCAI_ASSERT_EQUAL_ERROR( mNumDiagonals, mDlg.size() )
    // check column indexes in JA
    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
        ContextPtr loc = this->getContextPtr();
        validIndexes.getSupportedContext( loc );
        ReadAccess<IndexType> rJA( mJa, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_ERROR( validIndexes[ loc ]( rJA.get(), mNumValues, mNumColumns ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }
    // ToDo: check ILG[0] == mNumDiagonals, be careful about size of ILG
    // check descending values in ILG, DLG
    {
        static LAMAKernel<UtilKernelTrait::isSorted<IndexType> > isSorted;
        ContextPtr loc = this->getContextPtr();
        isSorted.getSupportedContext( loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        SCAI_CONTEXT_ACCESS( loc )
        bool ascending = false; // check for descending
        SCAI_ASSERT_ERROR( isSorted[ loc ]( rIlg.get(), mNumRows, ascending ),
                           *this << " @ " << msg << ": not descending values in ILG" )
        SCAI_ASSERT_ERROR( isSorted[ loc ]( rDlg.get(), mNumDiagonals, ascending ),
                           *this << " @ " << msg << ": not descending values in DLG" )
    }
    // both, ILG and DLG, must sum up to mNumValues
    {
        static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
        ContextPtr loc = this->getContextPtr();
        reduce.getSupportedContext( loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_EQUAL_ERROR( reduce[loc]( rIlg.get(), mNumRows, 0, utilskernel::binary::ADD ), mNumValues )
        SCAI_ASSERT_EQUAL_ERROR( reduce[loc]( rDlg.get(), mNumDiagonals, 0, utilskernel::binary::ADD ), mNumValues )
    }

    // check index values in Perm for out of range

    if ( mNumRows > 0 )
    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
        ContextPtr loc = this->getContextPtr();
        validIndexes.getSupportedContext( loc );
        ReadAccess<IndexType> rJA( mJa, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_ERROR( validIndexes[loc]( rPerm.get(), mNumRows, mNumRows ),
                           *this << " @ " << msg << ": illegel indexes in Perm" )
    }

    // check perm: no values out of range, but make sure that it is permutation, e.g. [ 0, 0] is illegal

    if ( mNumRows > 0 ) // very important as maxval would not work
    {
        ContextPtr loc = getContextPtr();
        // temporary array for inverse permutation, initialize with mNumRows
        HArray<IndexType> invPermArray( mNumRows, mNumRows );
        static LAMAKernel<UtilKernelTrait::setInversePerm> setInversePerm;
        static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
        ReadAccess<IndexType> rPerm( mPerm, loc );
        WriteAccess<IndexType> wInversePerm( invPermArray, loc );
        SCAI_CONTEXT_ACCESS( loc )
        // set inverse permutation, should overwrite all values 'mNumRows'
        setInversePerm[loc]( wInversePerm.get(), rPerm.get(), mNumRows );
        IndexType maxIndex = reduce[loc]( wInversePerm.get(), mNumRows, 0, utilskernel::binary::MAX );
        SCAI_ASSERT_ERROR( maxIndex < mNumRows, "Perm array does not cover all row indexes, #rows = " << mNumRows );
    }

    // Note: check is not exhaustive, e.g. it does not check for same column index in one row
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_INFO( logger, "set identity values with size = " << size )
    mNumRows = size;
    mNumColumns = size;
    mNumDiagonals = 1; // identity has exactly one diagonal
    mNumValues = mNumRows;
    ContextPtr prefLoc = this->getContextPtr();
    mValues.clear();  // invalidate all values
    mValues.resize( mNumValues );
    HArrayUtils::setScalar( mValues, ValueType( 1 ), utilskernel::binary::COPY, prefLoc );
    HArrayUtils::setOrder( mPerm, mNumRows, prefLoc );
    HArrayUtils::setOrder( mJa,  mNumRows, prefLoc );
    mDlg.clear();
    mDlg.resize( mNumDiagonals );
    HArrayUtils::setScalar( mDlg, mNumRows, utilskernel::binary::COPY, prefLoc );
    mIlg.clear();
    mIlg.resize( mNumRows );
    HArrayUtils::setScalar( mIlg, IndexType( 1 ), utilskernel::binary::COPY, prefLoc );
    mDiagonalProperty = true;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setupData( ContextPtr context )
{
    SCAI_LOG_INFO( logger, "setupData" )
    SCAI_ASSERT_EQUAL_ERROR( mIlg.size(), mNumRows )
    static LAMAKernel<UtilKernelTrait::getValue<IndexType> > getValue;
    static LAMAKernel<JDSKernelTrait::ilg2dlg> ilg2dlg;
    ContextPtr loc = context;
    getValue.getSupportedContext( loc, ilg2dlg );
    ReadAccess<IndexType> ilg( mIlg, loc );
    WriteOnlyAccess<IndexType> dlg( mDlg, loc, mNumDiagonals );
    mNumDiagonals = 0;
    mNumValues = 0;
    SCAI_CONTEXT_ACCESS( loc )

    if ( mNumRows )
    {
        mNumDiagonals = getValue[loc]( ilg.get(), 0 );
    }

    mNumValues = ilg2dlg[loc]( dlg.get(), mNumDiagonals, ilg.get(), mNumRows );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::sortRows( ContextPtr context )
{
    SCAI_LOG_INFO( logger, *this << "sortRows, #rows = " << mNumRows )
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<UtilKernelTrait::sort<IndexType> > sortRows;
    ContextPtr loc = context;
    reduce.getSupportedContext( loc, sortRows );
    // sort the rows according to the array ilg, take sorting over in perm
    WriteAccess<IndexType> ilg( mIlg, loc );
    WriteAccess<IndexType> perm( mPerm, loc );
    SCAI_CONTEXT_ACCESS( loc )
    // reduce with ABS_MAX returns 0 ( instead of -max ) for mNumRows == 0
    mNumDiagonals = reduce[loc]( ilg.get(), mNumRows, 0, utilskernel::binary::ABS_MAX );
    SCAI_LOG_INFO( logger, *this << "sortRows on " << *loc << ", #jagged diagonals = " << mNumDiagonals )
    bool descending = false;
    sortRows[loc]( perm.get(), ilg.get(), ilg.get(), mNumRows, descending );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<OtherValueType>* values,
    const ContextPtr context ) const
{
    SCAI_REGION( "Storage.JDS.buildCSR" )
    SCAI_LOG_INFO( logger,
                   "buildCSR<" << common::getScalarType<OtherValueType>() << ">"
                   << " from JDS<" << common::getScalarType<ValueType>() << ">" << " on " << *context )
    static LAMAKernel<UtilKernelTrait::setScatter<IndexType, IndexType> > setScatter;
    static LAMAKernel<JDSKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;
    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;
    static LAMAKernel<UtilKernelTrait::setInversePerm> setInversePerm;
    ContextPtr loc = context;
    setScatter.getSupportedContext( loc );
    getCSRValues.getSupportedContext( loc, sizes2offsets, setInversePerm );
    // now we are sure to have a loc where all kernel routines have been implemented
    ReadAccess<IndexType> rJdsPerm( mPerm, loc );
    ReadAccess<IndexType> rJdsILG( mIlg, loc );
    WriteOnlyAccess<IndexType> wCsrIA( ia, loc, mNumRows + 1 );
    SCAI_CONTEXT_ACCESS( loc )
    // rowValues[ perm[i] ] = ilg[i]
    setScatter[loc]( wCsrIA.get(), rJdsPerm.get(), true, rJdsILG.get(), utilskernel::binary::COPY, mNumRows );

    if ( ja == NULL || values == NULL )
    {
        wCsrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets[loc]( wCsrIA.get(), mNumRows );
    SCAI_ASSERT_EQUAL_DEBUG( numValues, mNumValues )
    SCAI_LOG_DEBUG( logger, "buildCSR from JDS with " << mNumValues << " values" )
    // temporary array for inverse permutation
    HArray<IndexType> invPermArray; // allows to find a CSR row in JDS rows
    WriteOnlyAccess<IndexType> wJdsInversePerm( invPermArray, loc, mNumRows );
    // compute the inverse permutation so that we find original row in JDS data
    setInversePerm[loc]( wJdsInversePerm.get(), rJdsPerm.get(), mNumRows );
    WriteOnlyAccess<IndexType> wCsrJA( *ja, loc, mNumValues );
    WriteOnlyAccess<OtherValueType> wCsrValues( *values, loc, mNumValues );
    ReadAccess<IndexType> rJdsDLG( mDlg, loc );
    ReadAccess<IndexType> rJdsJA( mJa, loc );
    ReadAccess<ValueType> rJdsValues( mValues, loc );
    // now we can convert JDS to CSR via interface
    getCSRValues[loc]( wCsrJA.get(), wCsrValues.get(), wCsrIA.get(), mNumRows, wJdsInversePerm.get(), rJdsILG.get(),
                       rJdsDLG.get(), rJdsJA.get(), rJdsValues.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr prefLoc )
{
    if ( ia.size() == numRows )
    {
        // offset array required
        HArray<IndexType> offsets;
        IndexType total = _MatrixStorage::sizes2offsets( offsets, ia, prefLoc );
        SCAI_ASSERT_EQUAL( numValues, total, "sizes do not sum to number of values" );
        setCSRDataImpl( numRows, numColumns, numValues, offsets, ja, values, prefLoc );
        return;
    }

    SCAI_REGION( "Storage.JDS.setCSR" )
    SCAI_LOG_INFO( logger,
                   "setCSRDataImpl<" << common::getScalarType<ValueType>() << "," << common::getScalarType<OtherValueType>() << ">" << ", shape is " << numRows << " x " << numColumns << ", #values for CSR = " << numValues )
    static LAMAKernel<CSRKernelTrait::offsets2sizes> offsets2sizes;
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    static LAMAKernel<JDSKernelTrait::setCSRValues<ValueType, OtherValueType> > setCSRValues;
    ContextPtr loc = this->getContextPtr();
    offsets2sizes.getSupportedContext( loc, setCSRValues );
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
        offsets2sizes[loc]( wIlg.get(), rCsrIA.get(), mNumRows );
        // set perm to identity permutation
        setOrder[loc]( wPerm.get(), mNumRows );
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
        setCSRValues[loc]( wJa.get(), wValues.get(), numRows, rPerm.get(), rIlg.get(), numDiagonals, rDlg.get(),
                           rCsrIA.get(), rCsrJA.get(), rCsrValues.get() );
    }
    this->resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setDIADataImpl(
    const IndexType /*numRows*/,
    const IndexType /*numColumns*/,
    const IndexType /*numDiagonals*/,
    const HArray<IndexType>& /*offsets*/,
    const HArray<OtherValueType>& /*values*/,
    const ContextPtr /*prefLoc*/ )
{
    COMMON_THROWEXCEPTION( "not yet implemeted" )
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

    if ( mNumRows > 0 )
    {
        // the arrays mIlg and mPerm need initialization
        static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;
        static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
        ContextPtr loc = this->getContextPtr();
        setOrder.getSupportedContext( loc, setVal );
        mNumRows = numRows;
        mNumColumns = numColumns;
        // we allocate at least ilg and perm with the correct value for a zero matrix
        SCAI_CONTEXT_ACCESS( loc )
        WriteOnlyAccess<IndexType> ilg( mIlg, loc, mNumRows );
        WriteOnlyAccess<IndexType> perm( mPerm, loc, mNumRows );
        setVal[loc]( ilg.get(), mNumRows, 0, utilskernel::binary::COPY );
        setOrder[loc]( perm.get(), mNumRows );
    }

    mDiagonalProperty = checkDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "JDSStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << mNumRows << " x " << mNumColumns
           << ", jd = " << mNumDiagonals << ", nnz = " << mNumValues << " )";
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<JDSKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );

    IndexType pos = getValuePos[loc]( i, j, mNumRows, dlg.get(), ilg.get(), perm.get(), ja.get() );

    ValueType val = 0;

    if ( pos != nIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumValues, "illegal value position for ( " << i << ", " << j << " )" );

        val = utilskernel::HArrayUtils::getVal<ValueType>( mValues, pos );
    }

    return val;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const utilskernel::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_DEBUG( logger, "set value (" << i << ", " << j << ")" )

    static LAMAKernel<JDSKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJa, loc );

    IndexType pos = getValuePos[loc]( i, j, mNumRows, dlg.get(), ilg.get(), perm.get(), ja.get() );

    if ( pos == nIndex )
    {
        COMMON_THROWEXCEPTION( "ELL storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumValues, "illegal value position for ( " << i << ", " << j << " )" );

    utilskernel::HArrayUtils::setValImpl( mValues, pos, val, op );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.JDS.timesVector" )
    SCAI_LOG_DEBUG( logger,
                    "Computing z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    if ( alpha == common::constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::binaryOpScalar1( result, beta, y, utilskernel::binary::MULT, this->getContextPtr() );
        return;
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )

    if ( beta != common::constants::ZERO )
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
    }

    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )
    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    SCAI_CONTEXT_ACCESS( loc )

    // this call will finish the computation, syncToken == NULL

    if ( beta != common::constants::ZERO )
    {
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, jdsPerm.get(), jdsILG.get(),
                         mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get() );
    }
    else
    {
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, NULL, mNumRows, jdsPerm.get(), jdsILG.get(),
                         mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get() );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.JDS.vectorTimesMatrix" )
    SCAI_LOG_DEBUG( logger,
                    "Computing z = " << alpha << " * x * A + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )
    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    ContextPtr loc = this->getContextPtr();

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumColumns );
        HArrayUtils::setScalar( result, ValueType( 0 ), utilskernel::binary::COPY, loc );
    }
    else
    {
        // Note: assignScaled will deal with
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
        HArrayUtils::binaryOpScalar1( result, beta, y, utilskernel::binary::MULT, loc );
    }

    // Step 2: result = alpha * x * this + 1 * result
    static LAMAKernel<JDSKernelTrait::normalGEVM<ValueType> > normalGEVM;
    normalGEVM.getSupportedContext( loc );
    SCAI_LOG_INFO( logger, *this << ": vectorTimesMatrix on " << *loc )
    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc, mNumColumns );
    SCAI_CONTEXT_ACCESS( loc )
    // this call will finish the computation, syncToken == NULL
    normalGEVM[loc]( wResult.get(), alpha, rX.get(), ValueType( 1 ), wResult.get(),
                     mNumColumns, jdsPerm.get(), jdsILG.get(),
                     mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
tasking::SyncToken* JDSStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    ContextPtr loc = getContextPtr();
    // For CUDA a solution using stream synchronization is more efficient than using a task
    SCAI_REGION( "Storage.JDS.timesVectorAsync" )
    SCAI_LOG_INFO( logger,
                   "Async start z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )
    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )
    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )
    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;
    normalGEMV.getSupportedContext( loc );
    common::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
    SCAI_ASYNCHRONOUS( syncToken.get() )
    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.
    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
    SCAI_CONTEXT_ACCESS( loc )
    // this call will only start the computation
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, jdsPerm.get(), jdsILG.get(),
                     mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get() );
    syncToken->pushRoutine( wResult.releaseDelayed() );
    syncToken->pushRoutine( jdsPerm.releaseDelayed() );
    syncToken->pushRoutine( jdsDLG.releaseDelayed() );
    syncToken->pushRoutine( jdsILG.releaseDelayed() );
    syncToken->pushRoutine( jdsJA.releaseDelayed() );
    syncToken->pushRoutine( jdsValues.releaseDelayed() );
    syncToken->pushRoutine( rX.releaseDelayed() );
    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
tasking::SyncToken* JDSStorage<ValueType>::vectorTimesMatrixAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    static LAMAKernel<JDSKernelTrait::normalGEVM<ValueType> > normalGEVM;
    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );
    // For CUDA a solution using stream synchronization is more efficient than using a task
    SCAI_REGION( "Storage.JDS.vectorTimesMatrixAsync" )
    SCAI_LOG_INFO( logger,
                   "Async start z = " << alpha << " * x * A + " << beta << " * y, with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )
    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumColumns )
    SCAI_LOG_INFO( logger, *this << ": matrixTimesVector on " << *loc )
    common::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
    SCAI_ASYNCHRONOUS( syncToken.get() )
    SCAI_CONTEXT_ACCESS( loc )
    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.
    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    if ( beta == scai::common::constants::ZERO )
    {
        // alias of result and y handled by correct order
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        // this call will only start the computation
        normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, NULL, mNumColumns, jdsPerm.get(), jdsILG.get(),
                         mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
    }
    else
    {
        // alias of result and y handled by correct order
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
        // this call will only start the computation
        normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumColumns, jdsPerm.get(), jdsILG.get(),
                         mNumDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
    }

    syncToken->pushRoutine( jdsPerm.releaseDelayed() );
    syncToken->pushRoutine( jdsDLG.releaseDelayed() );
    syncToken->pushRoutine( jdsILG.releaseDelayed() );
    syncToken->pushRoutine( jdsJA.releaseDelayed() );
    syncToken->pushRoutine( jdsValues.releaseDelayed() );
    syncToken->pushRoutine( rX.releaseDelayed() );
    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.JDS.jacobiIterate" )
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )
    static LAMAKernel<JDSKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL( mNumRows, mNumColumns, "storage must be square" )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, rhs.size() )
    // matrix must be square
    {
        ReadAccess<IndexType> jdsDlg( mDlg, loc );
        ReadAccess<IndexType> jdsIlg( mIlg, loc );
        ReadAccess<IndexType> jdsPerm( mPerm, loc );
        ReadAccess<IndexType> jdsJA( mJa, loc );
        ReadAccess<ValueType> jdsValues( mValues, loc );
        ReadAccess<ValueType> rOldSolution( oldSolution, loc );
        ReadAccess<ValueType> rRhs( rhs, loc );
        WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );
        SCAI_CONTEXT_ACCESS( loc )
        jacobi[loc]( wSolution.get(), mNumRows, jdsPerm.get(), jdsIlg.get(), mNumDiagonals, jdsDlg.get(), jdsJA.get(),
                     jdsValues.get(), rOldSolution.get(), rRhs.get(), omega );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
tasking::SyncToken* JDSStorage<ValueType>::jacobiIterateAsync(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.JDS.jacobiIterateAsync" )
    static LAMAKernel<JDSKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );

    if ( loc->getType() == Context::Host )
    {
        // On host we start directly a new task, avoids pushing accesses
        void ( JDSStorage::*jb )(
            HArray<ValueType>&,
            const HArray<ValueType>&,
            const HArray<ValueType>&,
            const ValueType omega ) const
        = &JDSStorage<ValueType>::jacobiIterate;
        using scai::common::bind;
        using scai::common::ref;
        using scai::common::cref;
        return new tasking::TaskSyncToken( bind( jb, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square
    common::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
    syncToken->setCurrent();
    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsJA( mJa, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    SCAI_CONTEXT_ACCESS( loc )
    jacobi[loc]( wSolution.get(), mNumRows, jdsPerm.get(), jdsILG.get(), mNumDiagonals, jdsDLG.get(), jdsJA.get(),
                 jdsValues.get(), rOldSolution.get(), rRhs.get(), omega );
    syncToken->pushRoutine( wSolution.releaseDelayed() );
    syncToken->pushRoutine( jdsDLG.releaseDelayed() );
    syncToken->pushRoutine( jdsILG.releaseDelayed() );
    syncToken->pushRoutine( jdsPerm.releaseDelayed() );
    syncToken->pushRoutine( jdsJA.releaseDelayed() );
    syncToken->pushRoutine( jdsValues.releaseDelayed() );
    syncToken->pushRoutine( rOldSolution.releaseDelayed() );
    syncToken->pushRoutine( rRhs.releaseDelayed() );
    syncToken->unsetCurrent();
    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const HArray<ValueType>& oldHaloSolution,
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
    const HArray<ValueType>* localDiagonal;
    shared_ptr<HArray<ValueType> > tmpLocalDiagonal;
    tmpLocalDiagonal = shared_ptr<HArray<ValueType> >( new HArray<ValueType>() );
    localStorage.getDiagonal( *tmpLocalDiagonal );
    localDiagonal = tmpLocalDiagonal.get();
    jacobiIterateHalo( localSolution, *localDiagonal, oldHaloSolution, omega );
}

/* ------------------------------------------------------------------------------------------------------------------ */
template<typename ValueType>
void JDSStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )
    SCAI_REGION( "Storage.JDS.jacobiIterateHalo" )
    static LAMAKernel<JDSKernelTrait::jacobiHalo<ValueType> > jacobiHalo;
    ContextPtr loc = this->getContextPtr();
    jacobiHalo.getSupportedContext( loc );
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
    jacobiHalo[loc]( wSolution.get(), mNumRows, diagonal.get(), mNumDiagonals, jdsHaloPerm.get(), jdsHaloIlg.get(),
                     jdsHaloDlg.get(), jdsHaloJA.get(), jdsHaloValues.get(), rOldHaloSolution.get(), omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType JDSStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )
    const IndexType n = mNumValues;

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    ContextPtr loc = this->getContextPtr();
    asum.getSupportedContext( loc );
    ReadAccess<ValueType> data( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    return asum[loc]( n, data.get(), 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType JDSStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )
    const IndexType n = mNumValues;

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = this->getContextPtr();
    dot.getSupportedContext( loc );
    ReadAccess<ValueType> data( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    return common::Math::sqrt( dot[loc]( n, data.get(), 1, data.get(), 1 ) );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )
    return mValues.maxNorm();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "JDSStorage of matrix " << mNumRows << " x " << mNumColumns;
    stream << ", #non-zero values = " << mNumValues << endl;
    ReadAccess<IndexType> perm( mPerm );
    ReadAccess<IndexType> ilg( mIlg );
    ReadAccess<IndexType> dlg( mDlg );
    ReadAccess<IndexType> ja( mJa );
    ReadAccess<ValueType> values( mValues );

    for ( IndexType ii = 0; ii < mNumRows; ii++ )
    {
        stream << "   row " << ii << " is original row " << perm[ii];
        stream << ", #non-zero values = " << ilg[ii] << endl;
        IndexType offset = ii;
        stream << "     column indexes = ";

        for ( IndexType d = 0; d < ilg[ii]; d++ )
        {
            stream << " " << ja[offset];
            offset += dlg[d];
        }

        stream << endl;
        offset = ii;
        stream << "     values   = ";

        for ( IndexType d = 0; d < ilg[ii]; d++ )
        {
            stream << " " << values[offset];
            offset += dlg[d];
        }

        stream << endl;
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::prefetch( const ContextPtr context ) const
{
    mDlg.prefetch( context );
    mIlg.prefetch( context );
    mPerm.prefetch( context );
    mJa.prefetch( context );
    mValues.prefetch( context );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const HArray<IndexType>& JDSStorage<ValueType>::getJA() const
{
    return mJa;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const HArray<IndexType>& JDSStorage<ValueType>::getPerm() const
{
    return mPerm;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const HArray<IndexType>& JDSStorage<ValueType>::getDlg() const
{
    return mDlg;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const HArray<IndexType>& JDSStorage<ValueType>::getIlg() const
{
    return mIlg;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const HArray<ValueType>& JDSStorage<ValueType>::getValues() const
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

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::swap( _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( getFormat(), other.getFormat(), "swap only for same storage format" )
    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(), "swap only for same value type" )

    // only in debug mode use the more expensive dynamic cast for verification

    SCAI_ASSERT_DEBUG( dynamic_cast<JDSStorage<ValueType>* >( &other ), "illegal storage to swap" )

    swapImpl( reinterpret_cast<JDSStorage<ValueType>& >( other ) );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::swapImpl( JDSStorage<ValueType>& other )
{
    MatrixStorage<ValueType>::swapMS( other ); // swap member variable of base class
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
    memoryUsage += sizeof( IndexType );
    memoryUsage += sizeof( IndexType );
    memoryUsage += sizeof( IndexType ) * mDlg.size();
    memoryUsage += sizeof( IndexType ) * mIlg.size();
    memoryUsage += sizeof( IndexType ) * mPerm.size();
    memoryUsage += sizeof( IndexType ) * mJa.size();
    memoryUsage += sizeof( ValueType ) * mValues.size();
    return memoryUsage;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>* JDSStorage<ValueType>::copy() const
{
    return new JDSStorage( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
JDSStorage<ValueType>* JDSStorage<ValueType>::newMatrixStorage() const
{
    common::unique_ptr<JDSStorage<ValueType> > storage( new JDSStorage<ValueType>() );
    storage->setContextPtr( this->getContextPtr() );
    return storage.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* JDSStorage<ValueType>::create()
{
    return new JDSStorage<ValueType>();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType JDSStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::JDS, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string JDSStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "JDSStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* JDSStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( JDSStorage, SCAI_NUMERIC_TYPES_HOST )

#define JDS_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                 \
    template void JDSStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                              \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                               \
    template void JDSStorage<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;              \
    template void JDSStorage<ValueType>::setRowImpl( const hmemo::HArray<OtherValueType>&, const IndexType,                \
            const utilskernel::binary::BinaryOp );                          \
    template void JDSStorage<ValueType>::getColumnImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;           \
    template void JDSStorage<ValueType>::setColumnImpl( const hmemo::HArray<OtherValueType>&, const IndexType,             \
            const utilskernel::binary::BinaryOp );                       \
    template void JDSStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                          \
    template void JDSStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                          \
    template void JDSStorage<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                                \
    template void JDSStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;                                               \
    template void JDSStorage<ValueType>::setDIADataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );

#define JDS_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, JDS_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( JDS_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef JDS_STORAGE_INST_LVL2
#undef JDS_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
