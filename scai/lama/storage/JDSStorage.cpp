/**
 * @file JDSStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/sparsekernel/JDSUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

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

using utilskernel::HArrayUtils;

using sparsekernel::JDSUtils;
using sparsekernel::CSRUtils;

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
    if ( other.getFormat() == Format::JDS )
    {
        // both storage have JDS format, use special method for it
        
        assignJDS( static_cast<const JDSStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRData( otherCSR.getNumRows(), otherCSR.getNumColumns(), 
                    otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues() );
    }
    else
    {
        HArray<IndexType>  csrIA;
        HArray<IndexType>  csrJA;
        HArray<ValueType>  csrValues;     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        // just a thought for optimization: use mIA, mJA, mValues instead of csrIA, csrJA, csrValues
        // but does not help much at all as resort of entries requires already temporaries.
        
        setCSRData( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues );
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

    JDSUtils::getRow( row, getNumColumns(), i,
                      mIlg, mDlg, mPerm, mJA, mValues, getContextPtr() ); 
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.JDS.getSparseRow" )

    SCAI_LOG_INFO( logger, "getSparseRow( " << i << " of " << getNumRows() << " )" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    HArray<IndexType> pos;  // positions in the array mJA, mValues

    JDSUtils::getRowPositions( pos, mIlg, mDlg, mPerm, i, getContextPtr() );

    // with entries in pos we can gather the column indexes and the values from jdsJA, jdsValues
 
    HArrayUtils::gather( jA, mJA, pos, BinaryOp::COPY, getContextPtr() );
    HArrayUtils::gather( values, mValues, pos, BinaryOp::COPY, getContextPtr() );

    SCAI_LOG_DEBUG( logger, "getSparseRow( " << i << " ) : jA = " << jA << ", values = " << values )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const
{   
    SCAI_REGION( "Storage.JDS.getSparseCol" )
    
    SCAI_LOG_INFO( logger, "getColumn( " << j << " ) of : " << *this )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )
    
    HArray<IndexType> pos;   // temparary array for the positions of the column entries

    JDSUtils::getColumnPositions( iA, pos, mIlg, mDlg, mPerm, mJA, j, getContextPtr() );

    HArrayUtils::gather( values, mValues, pos, BinaryOp::COPY, getContextPtr() );
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
    SCAI_ASSERT_EQ_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    JDSUtils::setRow( mValues, i, row, mIlg, mDlg, mPerm, mJA, op, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const BinaryOp op )
{
    SCAI_LOG_INFO( logger, "setColumn( " << j << " ) of : " << *this << " with column " << column )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), getNumRows(), "column array to small for set" )

    SCAI_REGION( "Storage.JDS.setCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    JDSUtils::getColumnPositions( rowIndexes, valuePos, mIlg, mDlg, mPerm, mJA, j, getContextPtr() );

    SCAI_LOG_INFO( logger, "setColumn( " << j << " ) updates " << rowIndexes.size() << " entries" )

    //  mValues[ pos ] op= column[row], scatter indexes into mValues are unique

    HArrayUtils::gather( colValues, column, rowIndexes, BinaryOp::COPY, getContextPtr() );
    HArrayUtils::scatter( mValues, valuePos, true, colValues, op, getContextPtr() );
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

    auto op = common::BinaryOp::MULT;

    JDSUtils::setRows( mValues, mIlg, mDlg, mPerm, diagonal, op, getContextPtr() );
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

    // check column indexes in JA

    SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mJA, getNumColumns(), getContextPtr() ),
                       this << " @ " << msg << ": illegal column indexes in JA"  );

    if  ( numRows > 0 )
    {
         IndexType maxValuesPerRow = mIlg[0];
         SCAI_ASSERT_EQ_ERROR( maxValuesPerRow, mDlg.size(), "serious inconsistency" )
    }

    // check descending values in ILG, DLG

    SCAI_ASSERT_ERROR( HArrayUtils::isSorted( mIlg, common::CompareOp::GE, getContextPtr() ),
                       *this << " @ " << msg << ": not decreasing values in ILG" )

    SCAI_ASSERT_ERROR( HArrayUtils::isSorted( mDlg, common::CompareOp::GE, getContextPtr() ),
                       *this << " @ " << msg << ": not decreasing values in DLG" )

    IndexType sumILG = HArrayUtils::reduce( mIlg, BinaryOp::ADD, getContextPtr() );
    IndexType sumDLG = HArrayUtils::reduce( mDlg, BinaryOp::ADD, getContextPtr() );

    SCAI_ASSERT_EQ_ERROR( sumILG, numValues, "inconsistent entries in jds ILG" )
    SCAI_ASSERT_EQ_ERROR( sumDLG, numValues, "inconsistent entries in jds DLG" )

    // check index values in Perm for out of range

    if ( numRows > 0 )
    {
        SCAI_ASSERT_ERROR( HArrayUtils::validIndexes( mPerm, getNumRows(), getContextPtr() ), "perm contains illegal indexes" )
    }

    // check perm: no values out of range, but make sure that it is permutation, e.g. [ 0, 0] is illegal

    if ( numRows > 0 ) // very important as maxval would not work
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
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void JDSStorage<ValueType>::buildCSRData(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    _HArray& csrValues ) const
{
    SCAI_REGION( "Storage.JDS.buildCSR" )

    SCAI_LOG_INFO( logger,
                   "buildCSR<" << csrValues.getValueType() << ">"
                   << " from JDS<" << common::getScalarType<ValueType>() << ">"
                   << " on " << *getContextPtr() << " ( preferred )" )

    if ( csrValues.getValueType() == this->getValueType() )
    {
        auto& castCSRValues = static_cast<HArray<ValueType>&>( csrValues );
        JDSUtils::convertJDS2CSR( csrIA, csrJA, castCSRValues, getNumRows(), getNumColumns(), 
                                  mIlg, mDlg, mPerm, mJA, mValues, getContextPtr() );
    }
    else
    {
        HArray<ValueType> tmpValues;  // use temporary for conversion of values
        JDSUtils::convertJDS2CSR( csrIA, csrJA, tmpValues, getNumRows(), getNumColumns(), 
                                  mIlg, mDlg, mPerm, mJA, mValues, getContextPtr() );
        HArrayUtils::_assign( csrValues, tmpValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::buildCSRSizes( hmemo::HArray<IndexType>& csrSizes ) const
{
   JDSUtils::buildRowSizes( csrSizes, mIlg, mPerm, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::setCSRData(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const _HArray& csrValues )
{
    SCAI_REGION( "Storage.ELL.setCSR" )

    if ( csrIA.size() == numRows )
    {
        HArray<IndexType> offsetIA;
        CSRUtils::sizes2offsets( offsetIA, csrIA, getContextPtr() );
        setCSRData( numRows, numColumns, offsetIA, csrJA, csrValues );
        return;
    }

    if ( csrValues.getValueType() != this->getValueType() )
    {
        SCAI_LOG_INFO( logger, "setCSRData<" << csrValues.getValueType() << ">, convert values to " << getValueType() )

        HArray<ValueType> sameTypeCSRValues;
        HArrayUtils::_assign( sameTypeCSRValues, csrValues, getContextPtr() );
        setCSRData( numRows, numColumns, csrIA, csrJA, sameTypeCSRValues );
        return;
    }

    SCAI_LOG_INFO( logger, "setCSRData<" << getValueType() << "> " << numRows << " x " << numColumns << ", nnz = " << csrJA.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    // jdsValues have same type, so we only have to cast

    const auto& sameTypeCSRValues = static_cast<const HArray<ValueType>&>( csrValues );

    JDSUtils::convertCSR2JDS( mIlg, mDlg, mPerm, mJA, mValues, numRows, numColumns, csrIA, csrJA, sameTypeCSRValues, getContextPtr() );
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

        HArrayUtils::setSameValue( mIlg, getNumRows(), IndexType( 0 ), getContextPtr() );
        HArrayUtils::setOrder( mPerm, getNumRows(), getContextPtr() );
    }
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

    IndexType pos = JDSUtils::getValuePos( i, j, mIlg, mDlg, mPerm, mJA, getContextPtr() );

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
void JDSStorage<ValueType>::setValue( 
    const IndexType i,
    const IndexType j,
    const ValueType val,
    const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_DEBUG( logger, "set value (" << i << ", " << j << ")" )

    IndexType pos = JDSUtils::getValuePos( i, j, mIlg, mDlg, mPerm, mJA, getContextPtr() );

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

    SCAI_LOG_INFO( logger,
                   "gemv ( op = " << op << ", async = " << async
                   << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    MatrixStorage<ValueType>::gemvCheck( alpha, x, beta, y, op );  // checks for correct sizes

    SyncToken* token = NULL;

    token = JDSUtils::gemv( result, alpha, x, beta, y,
                            getNumRows(), getNumColumns(), 
                            mIlg, mDlg, mPerm, mJA, mValues,
                            op, async, getContextPtr() );

    return token;
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

    bool async = false;

    JDSUtils::jacobi( solution, omega, oldSolution, rhs, 
                      mIlg, mDlg, mPerm, mJA, mValues, async, getContextPtr() );
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
 
    bool async = true;

    SyncToken* token = 

        JDSUtils::jacobi( solution, omega, oldSolution, rhs,
                          mIlg, mDlg, mPerm, mJA, mValues, async, getContextPtr() );

    return token;
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

    JDSUtils::jacobiHalo( localSolution, omega, oldHaloSolution, localDiagonal,
                          mIlg, mDlg, mPerm, mJA, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void JDSStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{   
    halo.halo2Global( mJA );
    _MatrixStorage::setDimension( getNumRows(), globalNumColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> JDSStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    const IndexType n = mValues.size();

    if ( n == 0 )
    {
        return RealType<ValueType>( 0 );
    }

    return HArrayUtils::l1Norm( mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> JDSStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    const IndexType n = mValues.size();

    if ( n == 0 )
    {
        return RealType<ValueType>( 0 );
    }

    return HArrayUtils::l2Norm( mValues, getContextPtr() );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
RealType<ValueType> JDSStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    return HArrayUtils::maxNorm( mValues, getContextPtr() );
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

} /* end namespace lama */

} /* end namespace scai */
