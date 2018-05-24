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
#include <scai/lama/storage/CSRStorage.hpp>

// local scai libraries
#include <scai/sparsekernel/JDSKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/JDSUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>
#include <functional>

using namespace scai::hmemo;

using std::shared_ptr;

namespace scai
{

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;

using sparsekernel::CSRKernelTrait;
using sparsekernel::JDSKernelTrait;
using sparsekernel::JDSUtils;

using tasking::SyncToken;

using common::BinaryOp;

namespace lama
{
// Allow for shared_ptr<ValueType> instead of std::shared_ptr<ValueType>


/* ------------------------------------------------------------------------------------------------------------------ */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, JDSStorage<ValueType>::logger, "MatrixStorage.JDSStorage" )

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( ContextPtr ctx ) : 

    MatrixStorage<ValueType>( 0, 0, ctx ),
    mDlg( ctx ),
    mIlg( ctx ),
    mPerm( ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "JDSStorage, matrix is 0 x 0." )

    _MatrixStorage::resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( IndexType numRows, IndexType numColumns, ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mDlg( ctx ),
    mIlg( numRows, IndexType( 0 ), ctx ),
    mPerm( ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    HArrayUtils::setOrder( mPerm, numRows, ctx );
    SCAI_LOG_DEBUG( logger, "JDSStorage<" << common::TypeTraits<ValueType>::id()
                            << ">( " << numRows << " x " << numColumns << " @ " << *ctx << " )" )

    _MatrixStorage::resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage(
    const IndexType numRows,
    const IndexType numColumns,
    HArray<IndexType> dlg,
    HArray<IndexType> ilg,
    HArray<IndexType> perm,
    HArray<IndexType> ja,
    HArray<ValueType> values,
    const ContextPtr ctx ) : 

    MatrixStorage<ValueType>( numRows, numColumns, ctx ), 
    mDlg( std::move( dlg ) ), 
    mIlg( std::move( ilg ) ), 
    mPerm( std::move( perm ) ), 
    mJA( std::move( ja ) ), 
    mValues( std::move( values ) )
{
    check( "JDSStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
    _MatrixStorage::resetDiagonalProperty();
    SCAI_LOG_INFO( logger, *this << ": constructed by JDS arrays dlg, ilg, .., values" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( JDSStorage<ValueType>&& other ) noexcept :

    MatrixStorage<ValueType>( std::move( other ) ),

    mDlg( std::move( other.mDlg ) ),
    mIlg( std::move( other.mIlg ) ),
    mPerm( std::move( other.mPerm ) ),
    mJA( std::move( other.mJA ) ),
    mValues( std::move( other.mValues ) )
{
    // no further checks as we assume a consistent and valid input JDS storage other
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setJDSData(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& dlg,
    const HArray<IndexType>& ilg,
    const HArray<IndexType>& perm,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    SCAI_ASSERT_EQUAL_ERROR( numRows, ilg.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows, perm.size() )
    SCAI_ASSERT_EQUAL_ERROR( values.size(), ja.size() )
    _MatrixStorage::setDimension( numRows, numColumns );
    ContextPtr loc = getContextPtr();
    HArrayUtils::assign( mDlg, dlg, loc );
    HArrayUtils::assign( mIlg, ilg, loc );
    HArrayUtils::assign( mPerm, perm, loc );
    HArrayUtils::assign( mJA, ja, loc );
    HArrayUtils::_assign( mValues, values, loc ); // supports type conversion
    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode
#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "JDSStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif
    this->resetDiagonalProperty();
    SCAI_LOG_INFO( logger, *this << ": set JDS by arrays dlg, ilg, .., values" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::JDSStorage( const JDSStorage<ValueType>& other ) : 

    MatrixStorage<ValueType>( other )
{
    assignJDS( other );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>& JDSStorage<ValueType>::operator=( const JDSStorage<ValueType>& other )
{
    assignJDS( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
JDSStorage<ValueType>& JDSStorage<ValueType>::operator=( JDSStorage<ValueType>&& other )
{
    // move of all member variables

    mDlg = std::move( other.mDlg );
    mIlg = std::move( other.mIlg );
    mPerm = std::move( other.mPerm );
    mJA = std::move( other.mJA );
    mValues = std::move( other.mValues );

    // call of move assignment for base class 

    MatrixStorage<ValueType>::moveImpl( std::move( other ) );

    return *this;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper
    
    mepr::StorageWrapper<JDSStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    ContextPtr ctx = getContextPtr();   // will force a valid copy in this context
    
    if ( other.getFormat() == Format::JDS )
    {
        // both storage have JDS format, use special method for it
        
        assignJDS( static_cast<const JDSStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRDataImpl( otherCSR.getNumRows(), otherCSR.getNumColumns(), 
                        otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues(), ctx );
    }
    else
    {
        HArray<IndexType>  csrIA( ctx );
        HArray<IndexType>  csrJA( ctx );
        HArray<ValueType>  csrValues( ctx );     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        // just a thought for optimization: use mIA, mJA, mValues instead of csrIA, csrJA, csrValues
        // but does not help much at all as resort of entries requires already temporaries.
        
        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues, ctx );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::assignJDS( const JDSStorage<OtherValueType>& other )
{
    if ( static_cast<const _MatrixStorage*>( &other ) == this )
    {
        SCAI_LOG_INFO( logger, typeName() << ": self assign, skipped, storage = " << other )
        return;
    }

    auto ctx = getContextPtr();

    // both storage have JDS format, we can just copy the corresponding arrays to the right context

    _MatrixStorage::_assign( other );     // assign member variables of base class

    HArrayUtils::assign( mIlg, other.getIlg(), ctx );
    HArrayUtils::assign( mDlg, other.getDlg(), ctx );
    HArrayUtils::assign( mPerm, other.getPerm(), ctx );
    HArrayUtils::assign( mJA, other.getJA(), ctx );
    HArrayUtils::assign( mValues, other.getValues(), ctx );

    _MatrixStorage::resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::splitUp(
    IndexType& numRows,
    IndexType& numColumns,
    hmemo::HArray<IndexType>& dlg,
    hmemo::HArray<IndexType>& ilg,
    hmemo::HArray<IndexType>& perm,
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values )
{
    // Take over all allocated data in the output arguments

    dlg    = std::move( mDlg );
    ilg    = std::move( mIlg );
    perm   = std::move( mPerm );
    ja     = std::move( mJA );
    values = std::move( mValues );

    _MatrixStorage::splitUp( numRows, numColumns );

    // this storage ends up in a valid zero storage
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::clear()
{
    _MatrixStorage::setDimension( 0, 0 );

    // clear all LAMA arrays used for this storage

    mDlg.clear();
    mIlg.clear();
    mPerm.clear();
    mJA.clear();
    mValues.clear();

    _MatrixStorage::resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
Format JDSStorage<ValueType>::getFormat() const
{
    return Format::JDS;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType JDSStorage<ValueType>::getNumValues() const
{
    return mValues.size();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType JDSStorage<ValueType>::getNumDiagonals() const
{
    return mDlg.size();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setDiagonal( const ValueType value )
{
    const IndexType numDiagonals = common::Math::min( getNumRows(), getNumColumns() );

    HArray<ValueType> diag( numDiagonals, value, getContextPtr() );

    setDiagonalV( diag );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    HArray<IndexType> diagonalPositions;

    IndexType numDiagonalsFound = JDSUtils::getDiagonalPositions(
        diagonalPositions, getNumRows(), getNumColumns(), mIlg, mDlg, mPerm, mJA, getContextPtr() );

    // as we have the number of found diagonals we have not to check for any invalidIndex

    SCAI_ASSERT_EQ_ERROR( diagonalPositions.size(), numDiagonalsFound,
                          "no diagonal property, some diagonal elements are missing" )

    SCAI_ASSERT_EQ_ERROR( diagonal.size(), diagonalPositions.size(), "diagonal has illegal size" )

    bool unique = true;  // there are no multiple diagonal entries

    HArrayUtils::scatter( mValues, diagonalPositions, unique, diagonal, common::BinaryOp::COPY, getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::getRow( HArray<ValueType>& row, const IndexType i ) const
{
    SCAI_REGION( "Storage.JDS.getRow" )

    SCAI_LOG_INFO( logger, "getRow with i = " << i )
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    static LAMAKernel<JDSKernelTrait::getRow<ValueType> > getRow;
    ContextPtr loc = this->getContextPtr();
    getRow.getSupportedContext( loc );
    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJA, loc );
    ReadAccess<ValueType> values( mValues, loc );
    WriteOnlyAccess<ValueType> wRow( row, loc, getNumColumns() );
    SCAI_CONTEXT_ACCESS( loc )
    getRow[loc]( wRow.get(), i, getNumColumns(), getNumRows(), perm.get(), ilg.get(), dlg.get(), ja.get(), values.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.JDS.getSparseRow" )

    SCAI_LOG_INFO( logger, "getSparseRow( " << i << " of " << getNumRows() << " )" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    static LAMAKernel<JDSKernelTrait::getValuePosRow> getValuePosRow;

    ContextPtr loc = this->getContextPtr();

    getValuePosRow.getSupportedContext( loc );

    HArray<IndexType> pos;  // positions in the array mJA, mValues

    {
        SCAI_CONTEXT_ACCESS( loc )

        // start with maximal possible size, is number of columns, resize later

        WriteOnlyAccess<IndexType> wPos( pos, loc, getNumColumns() );  

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );

        IndexType cnt = getValuePosRow[loc]( wPos.get(), i, getNumRows(),
                                             rIlg.get(), rDlg.get(), rPerm.get() );

        wPos.resize( cnt );
    }

    // with entries in pos we can gather the column indexes and the values from jdsJA, jdsValues
 
    HArrayUtils::gather( jA, mJA, pos, BinaryOp::COPY, loc );
    HArrayUtils::gather( values, mValues, pos, BinaryOp::COPY, loc );

    SCAI_LOG_DEBUG( logger, "getSparseRow( " << i << " ) : jA = " << jA << ", values = " << values )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const
{   
    SCAI_REGION( "Storage.JDS.getSparseCol" )
    
    SCAI_LOG_INFO( logger, "getColumn( " << j << " ) of : " << *this )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )
    
    static LAMAKernel<JDSKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( iA, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, getNumRows() );

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rJa( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j, getNumRows(),
                                             rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    HArrayUtils::gather( values, mValues, valuePos, BinaryOp::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.JDS.getDenseCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, getNumRows(), colValues, rowIndexes, ValueType ( 0 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const BinaryOp op )
{
    SCAI_REGION( "Storage.JDS.setRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    SCAI_LOG_INFO( logger, "setRowImpl( i = " << i << " )" )

    static LAMAKernel<JDSKernelTrait::setRow<ValueType> > setRow;

    ContextPtr loc = this->getContextPtr();
    setRow.getSupportedContext( loc );

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJA, loc );
    WriteAccess<ValueType> values( mValues, loc );
    ReadAccess<ValueType> rRow( row, loc );
    SCAI_CONTEXT_ACCESS( loc )
    setRow[loc]( values.get(), i, getNumColumns(), getNumRows(), perm.get(), ilg.get(), dlg.get(), ja.get(), rRow.get(), op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const BinaryOp op )
{
    SCAI_LOG_INFO( logger, "setColumn( " << j << " ) of : " << *this << " with column " << column )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), getNumRows(), "column array to small for set" )

    SCAI_REGION( "Storage.JDS.setCol" )

    static LAMAKernel<JDSKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, getNumRows() );

        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rJa( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j, getNumRows(),
                                             rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    SCAI_LOG_INFO( logger, "setColumn( " << j << " ) updates " << rowIndexes.size() << " entries" )

    //  mValues[ pos ] op= column[row]

    HArrayUtils::gather( colValues, column, rowIndexes, BinaryOp::COPY, loc );
    HArrayUtils::scatter( mValues, valuePos, true, colValues, op, loc );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    JDSUtils::getDiagonal( diagonal, getNumRows(), getNumColumns(), mIlg, mDlg, mPerm, mJA, mValues, getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::scale( const ValueType value )
{
    SCAI_LOG_INFO( logger, "scale with value = " << value )
    HArrayUtils::compute( mValues, mValues, BinaryOp::MULT, value, this->getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */


template<typename ValueType>
void JDSStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::UnaryOp::CONJ, this->getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::scaleRows( const HArray<ValueType>& diagonal )
{
    SCAI_LOG_INFO( logger, "scaleRows" )
    static LAMAKernel<JDSKernelTrait::scaleRows<ValueType> > jdsScaleRows;
    ContextPtr loc = this->getContextPtr();
    jdsScaleRows.getSupportedContext( loc );
    ReadAccess<ValueType> rDiagonal( diagonal, loc );
    ReadAccess<IndexType> rPerm( mPerm, loc );
    ReadAccess<IndexType> rIlg( mIlg, loc );
    ReadAccess<IndexType> rDlg( mDlg, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    jdsScaleRows[loc]( wValues.get(), getNumRows(), rPerm.get(), rIlg.get(), rDlg.get(), rDiagonal.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
bool JDSStorage<ValueType>::checkDiagonalProperty() const
{
    SCAI_LOG_DEBUG( logger, "checkDiagonalProperty: " << *this )

    HArray<IndexType> diagonalPositions;

    IndexType numDiagonalsFound = JDSUtils::getDiagonalPositions(
        diagonalPositions, getNumRows(), getNumColumns(), mIlg, mDlg, mPerm, mJA, getContextPtr() );

    return numDiagonalsFound == diagonalPositions.size();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::check( const char* msg ) const
{
    SCAI_LOG_DEBUG( logger, "check at " << *getContextPtr() << ", msg = " << msg )

    IndexType numRows = getNumRows();

    SCAI_ASSERT_EQUAL_ERROR( numRows, mIlg.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows, mPerm.size() )

    SCAI_ASSERT_EQUAL_ERROR( mValues.size(), mJA.size() )

    IndexType numValues    = mValues.size();
    IndexType numDiagonals = mDlg.size();

    // check column indexes in JA
    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
        ContextPtr loc = this->getContextPtr();
        validIndexes.getSupportedContext( loc );
        ReadAccess<IndexType> rJA( mJA, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_ERROR( validIndexes[ loc ]( rJA.get(), numValues, getNumColumns() ),
                           *this << " @ " << msg << ": illegel indexes in JA" )
    }
    // ToDo: check ILG[0] == numDiagonals, be careful about size of ILG
    // check descending values in ILG, DLG
    {
        static LAMAKernel<UtilKernelTrait::isSorted<IndexType> > isSorted;
        ContextPtr loc = this->getContextPtr();
        isSorted.getSupportedContext( loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_ERROR( isSorted[ loc ]( rIlg.get(), numRows, common::CompareOp::GE ),
                           *this << " @ " << msg << ": not decreasing values in ILG" )
        SCAI_ASSERT_ERROR( isSorted[ loc ]( rDlg.get(), numDiagonals, common::CompareOp::GE ),
                           *this << " @ " << msg << ": not decreasing values in DLG" )
    }
    // both, ILG and DLG, must sum up to num values
    {
        static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
        ContextPtr loc = this->getContextPtr();
        reduce.getSupportedContext( loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_EQUAL_ERROR( reduce[loc]( rIlg.get(), getNumRows(), 0, BinaryOp::ADD ), mValues.size() )
        SCAI_ASSERT_EQUAL_ERROR( reduce[loc]( rDlg.get(), numDiagonals, 0, BinaryOp::ADD ), mValues.size() )
    }

    // check index values in Perm for out of range

    if ( getNumRows() > 0 )
    {
        static LAMAKernel<UtilKernelTrait::validIndexes> validIndexes;
        ContextPtr loc = this->getContextPtr();
        validIndexes.getSupportedContext( loc );
        ReadAccess<IndexType> rJA( mJA, loc );
        ReadAccess<IndexType> rPerm( mPerm, loc );
        SCAI_CONTEXT_ACCESS( loc )
        SCAI_ASSERT_ERROR( validIndexes[loc]( rPerm.get(), getNumRows(), getNumRows() ),
                           *this << " @ " << msg << ": illegel indexes in Perm" )
    }

    // check perm: no values out of range, but make sure that it is permutation, e.g. [ 0, 0] is illegal

    if ( getNumRows() > 0 ) // very important as maxval would not work
    {
        ContextPtr loc = getContextPtr();

        HArray<IndexType> invPermArray;

        // this operation fails if perm is invalid

        HArrayUtils::inversePerm( invPermArray, mPerm );
    }

    // Note: check is not exhaustive, e.g. it does not check for same column index in one row
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_INFO( logger, "set identity values with size = " << size )

    _MatrixStorage::setDimension( size, size );

    ContextPtr ctx = this->getContextPtr();

    HArrayUtils::setSameValue( mValues, size, ValueType( 1 ), ctx );
    HArrayUtils::setOrder( mPerm, size, ctx );
    HArrayUtils::setOrder( mJA,  size, ctx );
    HArrayUtils::setSameValue( mIlg, size, IndexType( 1 ), ctx );

    mDlg = HArray<IndexType>( { getNumRows() }, ctx );  // one diagonal

    mDiagonalProperty = true;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    SCAI_LOG_INFO( logger, "set identity values with size = " << size )

    _MatrixStorage::setDimension( size, size );

    ContextPtr ctx = this->getContextPtr();

    HArrayUtils::setArray( mValues, diagonal, common::BinaryOp::COPY, ctx );
    HArrayUtils::setOrder( mPerm, size, ctx );
    HArrayUtils::setOrder( mJA,  size, ctx );
    HArrayUtils::setSameValue( mIlg, size, IndexType( 1 ), ctx );

    mDlg = HArray<IndexType>( { getNumRows() }, ctx );  // one diagonal

    mDiagonalProperty = true;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
IndexType JDSStorage<ValueType>::setupDiagonals()
{
    IndexType numRows = getNumRows();

    SCAI_ASSERT_EQUAL_ERROR( mIlg.size(), numRows )

    SCAI_LOG_INFO( logger, "setupDiagonals" )

    if ( 0 == numRows )
    {
        mDlg.clear();
        return 0;
    }

    IndexType numDiagonals = mIlg[0];  // is maximal number of entries in row

    static LAMAKernel<JDSKernelTrait::ilg2dlg> ilg2dlg;

    ContextPtr loc = getContextPtr();

    ilg2dlg.getSupportedContext( loc, ilg2dlg );

    ReadAccess<IndexType> ilg( mIlg, loc );
    WriteOnlyAccess<IndexType> dlg( mDlg, loc, numDiagonals );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType numValues = ilg2dlg[loc]( dlg.get(), numDiagonals, ilg.get(), getNumRows() );

    return numValues;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::sortRows()
{
    SCAI_LOG_INFO( logger, *this << "sortRows, #rows = " << getNumRows() )
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<UtilKernelTrait::sort<IndexType> > sortRowsKernel;
    ContextPtr loc = getContextPtr();
    reduce.getSupportedContext( loc, sortRowsKernel );
    // sort the rows according to the array ilg, take sorting over in perm
    WriteAccess<IndexType> ilg( mIlg, loc );
    WriteAccess<IndexType> perm( mPerm, loc );
    SCAI_CONTEXT_ACCESS( loc )
    // reduce with ABS_MAX returns 0 ( instead of -max ) for getNumRows() == 0

    SCAI_LOG_INFO( logger, *this << "sortRows on " << *loc )

    bool descending = false;
    sortRowsKernel[loc]( perm.get(), ilg.get(), ilg.get(), getNumRows(), descending );

    // Note: first entry is maximal, contains number of diagonals
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
    WriteOnlyAccess<IndexType> wCsrIA( ia, loc, getNumRows() + 1 );
    SCAI_CONTEXT_ACCESS( loc )
    // rowValues[ perm[i] ] = ilg[i]
    setScatter[loc]( wCsrIA.get(), rJdsPerm.get(), true, rJdsILG.get(), BinaryOp::COPY, getNumRows() );

    if ( ja == NULL || values == NULL )
    {
        wCsrIA.resize( getNumRows() );
        return;
    }

    IndexType numValues = sizes2offsets[loc]( wCsrIA.get(), getNumRows() );

    SCAI_ASSERT_EQ_DEBUG( numValues, mValues.size(), "row sizes do not sum up to number of nnz entries" )

    SCAI_LOG_DEBUG( logger, "buildCSR from JDS with " << numValues << " values" )
    // temporary array for inverse permutation
    HArray<IndexType> invPermArray; // allows to find a CSR row in JDS rows
    WriteOnlyAccess<IndexType> wJdsInversePerm( invPermArray, loc, getNumRows() );
    // compute the inverse permutation so that we find original row in JDS data
    setInversePerm[loc]( wJdsInversePerm.get(), rJdsPerm.get(), getNumRows() );
    WriteOnlyAccess<IndexType> wCsrJA( *ja, loc, numValues );
    WriteOnlyAccess<OtherValueType> wCsrValues( *values, loc, numValues );
    ReadAccess<IndexType> rJdsDLG( mDlg, loc );
    ReadAccess<IndexType> rJdsJA( mJA, loc );
    ReadAccess<ValueType> rJdsValues( mValues, loc );
    // now we can convert JDS to CSR via interface
    getCSRValues[loc]( wCsrJA.get(), wCsrValues.get(), wCsrIA.get(), getNumRows(), wJdsInversePerm.get(), rJdsILG.get(),
                       rJdsDLG.get(), rJdsJA.get(), rJdsValues.get() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void JDSStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr prefLoc )
{
    IndexType numValues = ja.size();

    if ( ia.size() == numRows )
    {
        // offset array required
        HArray<IndexType> offsets;
        IndexType total = _MatrixStorage::sizes2offsets( offsets, ia, prefLoc );
        SCAI_ASSERT_EQUAL( numValues, total, "sizes do not sum to number of values" );
        setCSRDataImpl( numRows, numColumns, offsets, ja, values, prefLoc );
        return;
    }

    SCAI_REGION( "Storage.JDS.setCSR" )
    SCAI_LOG_INFO( logger,
                   "setCSRDataImpl<" << common::getScalarType<ValueType>() << "," << common::getScalarType<OtherValueType>() << ">" << ", shape is " << numRows << " x " << numColumns << ", #values for CSR = " << ja.size() )
    static LAMAKernel<CSRKernelTrait::offsets2sizes> offsets2sizes;
    static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
    static LAMAKernel<JDSKernelTrait::setCSRValues<ValueType, OtherValueType> > setCSRValues;
    ContextPtr loc = this->getContextPtr();
    offsets2sizes.getSupportedContext( loc, setCSRValues );
    ReadAccess<IndexType> rCsrIA( ia, loc );
    ReadAccess<IndexType> rCsrJA( ja, loc );
    ReadAccess<OtherValueType> rCsrValues( values, loc );
    _MatrixStorage::setDimension( numRows, numColumns );
    // Step 1: fill up the array ilg and perm, detect diagonal property in csr data
    mDiagonalProperty = true; // will be set to false if not valid in one row
    {
        WriteOnlyAccess<IndexType> wIlg( mIlg, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wPerm( mPerm, loc, getNumRows() );
        SCAI_CONTEXT_ACCESS( loc )
        // ilg willl contain the sizes of each row
        offsets2sizes[loc]( wIlg.get(), rCsrIA.get(), getNumRows() );
        // set perm to identity permutation
        setOrder[loc]( wPerm.get(), getNumRows() );
    }

    sortRows(); // sorts ilg and builds perm

    numValues = setupDiagonals(); // sets dlg

    SCAI_ASSERT_EQ_ERROR( numValues, ja.size(), "sum of row sizes does not match size of ja" )
    SCAI_ASSERT_EQ_ERROR( numValues, values.size(), "sum of row sizes does not match size of values." )
 
    IndexType numDiagonals = mDlg.size();

    {
        ReadAccess<IndexType> rPerm( mPerm, loc );
        ReadAccess<IndexType> rIlg( mIlg, loc );
        ReadAccess<IndexType> rDlg( mDlg, loc );
        WriteOnlyAccess<ValueType> wValues( mValues, loc, numValues );
        WriteOnlyAccess<IndexType> wJa( mJA, loc, numValues );
        SCAI_CONTEXT_ACCESS( loc )
        setCSRValues[loc]( wJa.get(), wValues.get(), numRows, rPerm.get(), rIlg.get(), numDiagonals, rDlg.get(),
                           rCsrIA.get(), rCsrJA.get(), rCsrValues.get() );
    }

    _MatrixStorage::resetDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
JDSStorage<ValueType>::~JDSStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~JDSStorage( " << getNumRows() << " x " << getNumColumns() <<  " )" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::purge()
{
    _MatrixStorage::setDimension( 0, 0 );

    mDiagonalProperty = checkDiagonalProperty();
    mDlg.purge();
    mIlg.purge();
    mPerm.purge();
    mJA.purge();
    mValues.purge();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate JDS sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    _MatrixStorage::setDimension( numRows, numColumns );

    if ( getNumRows() > 0 )
    {
        // the arrays mIlg and mPerm need initialization
        static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;
        static LAMAKernel<UtilKernelTrait::setOrder<IndexType> > setOrder;
        ContextPtr loc = this->getContextPtr();
        setOrder.getSupportedContext( loc, setVal );

        // we allocate at least ilg and perm with the correct value for a zero matrix
        SCAI_CONTEXT_ACCESS( loc )
        WriteOnlyAccess<IndexType> ilg( mIlg, loc, getNumRows() );
        WriteOnlyAccess<IndexType> perm( mPerm, loc, getNumRows() );
        setVal[loc]( ilg.get(), getNumRows(), 0, BinaryOp::COPY );
        setOrder[loc]( perm.get(), getNumRows() );
    }

    mDiagonalProperty = checkDiagonalProperty();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "JDSStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << getNumRows() << " x " << getNumColumns()
           << ", jd = " << mDlg.size() << ", nnz = " << mValues.size() << " )";
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ValueType JDSStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<JDSKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, getNumRows(), ilg.get(), dlg.get(), perm.get(), ja.get() );

    ValueType val = 0;

    if ( pos != invalidIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mValues.size(), "illegal value position for ( " << i << ", " << j << " )" );

        val = utilskernel::HArrayUtils::getVal<ValueType>( mValues, pos );
    }

    return val;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_DEBUG( logger, "set value (" << i << ", " << j << ")" )

    static LAMAKernel<JDSKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> dlg( mDlg, loc );
    ReadAccess<IndexType> ilg( mIlg, loc );
    ReadAccess<IndexType> perm( mPerm, loc );
    ReadAccess<IndexType> ja( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, getNumRows(), ilg.get(), dlg.get(), perm.get(), ja.get() );

    if ( pos == invalidIndex )
    {
        COMMON_THROWEXCEPTION( "ELL storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mValues.size(), "illegal value position for ( " << i << ", " << j << " )" );

    utilskernel::HArrayUtils::setVal( mValues, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* JDSStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_REGION( "Storage.JDS.gemv" )

    const IndexType nSource = common::isTranspose( op ) ? getNumRows() : getNumColumns();
    const IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    IndexType numDiagonals = mDlg.size();

    SCAI_LOG_INFO( logger,
                   "GEMV ( op = " << op << ", async = " << async
                   << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::Constants::ZERO || ( numDiagonals == 0 ) )
    {
        // so we just have result = beta * y, will be done synchronously

        if ( beta == common::Constants::ZERO )
        {
            HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), this->getContextPtr() );
        }
        else
        {
            HArrayUtils::compute( result, beta, BinaryOp::MULT, y, this->getContextPtr() );
        }

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

    SCAI_ASSERT_EQUAL_ERROR( x.size(), nSource )
    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsJA( mJA, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    if ( beta != common::Constants::ZERO )
    {
        SCAI_ASSERT_EQ_ERROR( y.size(), nTarget, "y has illegal size" )

        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );  // result might be aliased to y

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), 
                         getNumRows(), getNumColumns(), jdsPerm.get(), jdsILG.get(),
                         numDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), op );
        if ( async )
        {
            syncToken->pushRoutine( rY.releaseDelayed() );
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }
    else
    {
        // do not access y at all

        WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL,
                         getNumRows(), getNumColumns(), jdsPerm.get(), jdsILG.get(),
                         numDiagonals, jdsDLG.get(), jdsJA.get(), jdsValues.get(), op );

        if ( async )
        {
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }

    if ( async )
    {
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( jdsPerm.releaseDelayed() );
        syncToken->pushRoutine( jdsDLG.releaseDelayed() );
        syncToken->pushRoutine( jdsILG.releaseDelayed() );
        syncToken->pushRoutine( jdsJA.releaseDelayed() );
        syncToken->pushRoutine( jdsValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    bool async = false; // synchronously execution, no SyncToken required
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
tasking::SyncToken* JDSStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y, 
    const common::MatrixOp op ) const
{
    bool async = true;
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token, "NULL token not allowed for asynchronous execution gemv, alpha = " << alpha << ", beta = " << beta )
    return token;
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

    SCAI_ASSERT_EQUAL( getNumRows(), getNumColumns(), "storage must be square" )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), rhs.size() )

    IndexType numDiagonals = mDlg.size();

    // matrix must be square
    {
        ReadAccess<IndexType> jdsDlg( mDlg, loc );
        ReadAccess<IndexType> jdsIlg( mIlg, loc );
        ReadAccess<IndexType> jdsPerm( mPerm, loc );
        ReadAccess<IndexType> jdsJA( mJA, loc );
        ReadAccess<ValueType> jdsValues( mValues, loc );
        ReadAccess<ValueType> rOldSolution( oldSolution, loc );
        ReadAccess<ValueType> rRhs( rhs, loc );
        WriteOnlyAccess<ValueType> wSolution( solution, loc, getNumRows() );
        SCAI_CONTEXT_ACCESS( loc )
        jacobi[loc]( wSolution.get(), getNumRows(), jdsPerm.get(), jdsIlg.get(), numDiagonals, jdsDlg.get(), jdsJA.get(),
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

    if ( loc->getType() == common::ContextType::Host )
    {
        // On host we start directly a new task, avoids pushing accesses
        void ( JDSStorage::*jb )(
            HArray<ValueType>&,
            const HArray<ValueType>&,
            const HArray<ValueType>&,
            const ValueType omega ) const
        = &JDSStorage<ValueType>::jacobiIterate;
        using std::bind;
        using std::ref;
        using std::cref;
        return new tasking::TaskSyncToken( bind( jb, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
    }

    // For CUDA a solution using stream synchronization is more efficient than using a task
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), getNumColumns() )

    IndexType numDiagonals = mDlg.size();

    // matrix must be square
    std::unique_ptr<tasking::SyncToken> syncToken( loc->getSyncToken() );
    syncToken->setCurrent();
    WriteOnlyAccess<ValueType> wSolution( solution, loc, getNumRows() );
    ReadAccess<IndexType> jdsDLG( mDlg, loc );
    ReadAccess<IndexType> jdsILG( mIlg, loc );
    ReadAccess<IndexType> jdsPerm( mPerm, loc );
    ReadAccess<IndexType> jdsJA( mJA, loc );
    ReadAccess<ValueType> jdsValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    SCAI_CONTEXT_ACCESS( loc )
    jacobi[loc]( wSolution.get(), getNumRows(), jdsPerm.get(), jdsILG.get(), numDiagonals, jdsDLG.get(), jdsJA.get(),
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
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for halo matrix data." )
    SCAI_REGION( "Storage.JDS.jacobiIterateHalo" )
    static LAMAKernel<JDSKernelTrait::jacobiHalo<ValueType> > jacobiHalo;
    ContextPtr loc = this->getContextPtr();
    jacobiHalo.getSupportedContext( loc );
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumColumns(), oldHaloSolution.size() )
    WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
    ReadAccess<ValueType> diagonal( localDiagonal, loc );
    ReadAccess<IndexType> jdsHaloPerm( mPerm, loc );
    ReadAccess<IndexType> jdsHaloIlg( mIlg, loc );
    ReadAccess<IndexType> jdsHaloDlg( mDlg, loc );
    ReadAccess<IndexType> jdsHaloJA( mJA, loc );
    ReadAccess<ValueType> jdsHaloValues( mValues, loc );
    ReadAccess<ValueType> rOldHaloSolution( oldHaloSolution, loc );

    IndexType numDiagonals = mDlg.size();

    SCAI_CONTEXT_ACCESS( loc )
    jacobiHalo[loc]( wSolution.get(), getNumRows(), diagonal.get(), numDiagonals, jdsHaloPerm.get(), jdsHaloIlg.get(),
                     jdsHaloDlg.get(), jdsHaloJA.get(), jdsHaloValues.get(), rOldHaloSolution.get(), omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{   
    halo.halo2Global( mJA );
    _MatrixStorage::setDimension( getNumRows(), globalNumColumns );
    _MatrixStorage::resetDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> JDSStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    const IndexType n = mValues.size();

    if ( n == 0 )
    {
        return static_cast<ValueType>( 0 );
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
RealType<ValueType> JDSStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    const IndexType n = mValues.size();

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
RealType<ValueType> JDSStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )
    return HArrayUtils::maxNorm( mValues );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "JDSStorage of matrix " << getNumRows() << " x " << getNumColumns();
    stream << ", #non-zero values = " << mValues.size() << endl;
    ReadAccess<IndexType> perm( mPerm );
    ReadAccess<IndexType> ilg( mIlg );
    ReadAccess<IndexType> dlg( mDlg );
    ReadAccess<IndexType> ja( mJA );
    ReadAccess<ValueType> values( mValues );

    for ( IndexType ii = 0; ii < getNumRows(); ii++ )
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
    mJA.prefetch( context );
    mValues.prefetch( context );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
const HArray<IndexType>& JDSStorage<ValueType>::getJA() const
{
    return mJA;
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
    mJA.wait();
    mValues.wait();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::swap( JDSStorage<ValueType>& other )
{
    // swap base class

    MatrixStorage<ValueType>::swap( other ); // swap member variable of base class

    // swap my member variables

    mDlg.swap( other.mDlg );
    mIlg.swap( other.mIlg );
    mPerm.swap( other.mPerm );
    mJA.swap( other.mJA );
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
    memoryUsage += sizeof( IndexType ) * mJA.size();
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
JDSStorage<ValueType>* JDSStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    std::unique_ptr<JDSStorage<ValueType> > storage( new JDSStorage<ValueType>( getContextPtr() ) );
    storage->allocate( numRows, numColumns );
    return storage.release();
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

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

template<typename ValueType>
const char* JDSStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType JDSStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::JDS, common::getScalarType<ValueType>() );
}

template<typename ValueType>
MatrixStorageCreateKeyType JDSStorage<ValueType>::getCreateValue() const
{
    return createValue();
}

template<typename ValueType>
_MatrixStorage* JDSStorage<ValueType>::create()
{
    return new JDSStorage<ValueType>();
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( JDSStorage, SCAI_NUMERIC_TYPES_HOST )

#define JDS_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                 \
    template void JDSStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType,                                 \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                              \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                               \
    template void JDSStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;

#define JDS_STORAGE_INST_LVL1( ValueType )                                                                                 \
    SCAI_COMMON_LOOP_LVL2( ValueType, JDS_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( JDS_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef JDS_STORAGE_INST_LVL2
#undef JDS_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
