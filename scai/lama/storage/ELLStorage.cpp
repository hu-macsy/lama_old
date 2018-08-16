/**
 * @file ELLStorage.cpp
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
 * @brief Instantitions for template class ELLStorage.
 * @author Lauretta Schubert
 * @date 25.05.2011
 */

// hpp
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/ELLKernelTrait.hpp>
#include <scai/sparsekernel/ELLUtils.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/hmemo.hpp>

#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tasking/NoSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>
#include <functional>

using std::unique_ptr;
using std::shared_ptr;

namespace scai
{

using tasking::SyncToken;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;

using sparsekernel::ELLKernelTrait;
using sparsekernel::ELLUtils;
using sparsekernel::CSRUtils;

using namespace tasking;
using namespace hmemo;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, ELLStorage<ValueType>::logger, "MatrixStorage.ELLStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( ContextPtr ctx ) : 

    MatrixStorage<ValueType>( 0, 0, ctx ), 
    mNumValuesPerRow( 0 ),
    mIA( ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "ELLStorage, default constructor for zero matrix." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( IndexType numRows, IndexType numColumns, ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mNumValuesPerRow( 0 ),
    mIA( numRows, IndexType( 0 ), ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "COOStorage for matrix " << getNumRows()
                             << " x " << getNumColumns() << ", no non-zero elements @ " << *ctx )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    HArray<IndexType> ia,
    HArray<IndexType> ja,
    HArray<ValueType> values,
    ContextPtr ctx ) : 

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mNumValuesPerRow( numValuesPerRow )
{
    // make some checks before the data is moved

    SCAI_ASSERT_EQUAL_ERROR( numRows, ia.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows * numValuesPerRow, ja.size() )
    SCAI_ASSERT_EQUAL_ERROR( numRows * numValuesPerRow, values.size() )

    mIA = std::move( ia );
    mJA = std::move( ja );
    mValues = std::move( values );

    fillValues();

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode
#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "ELLStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif

    SCAI_LOG_INFO( logger, *this << ": set ELLPACK by arrays ia, ja, values" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( const ELLStorage<ValueType>& other ) : 

    MatrixStorage<ValueType>( other )

{
    SCAI_LOG_INFO( logger, "copy constructor # other = " << other )

    // call the assignment operator

    assignELL( other );
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( ELLStorage<ValueType>&& other ) noexcept :

    MatrixStorage<ValueType>( std::move( other ) ),

    mNumValuesPerRow( other.mNumValuesPerRow ), 

    mIA( std::move( other.mIA ) ),
    mJA( std::move( other.mJA ) ),
    mValues( std::move( other.mValues ) )

{
    // no further checks as we assume other to be a consistent and valid input ELL storage 

    other.mNumValuesPerRow = 0;   // just to make sure that there remains a valid zero storage

    SCAI_LOG_DEBUG( logger, "moved ELLStorage, this = " << *this << ", other = " << other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>& ELLStorage<ValueType>::operator=( const ELLStorage<ValueType>& other )
{
    assignELL( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>& ELLStorage<ValueType>::operator=( ELLStorage<ValueType>&& other )
{
    // move of all member variables

    mIA = std::move( other.mIA );
    mJA = std::move( other.mJA );
    mValues = std::move( other.mValues );

    mNumValuesPerRow = other.mNumValuesPerRow;

    // call of move assignment for base class 

    MatrixStorage<ValueType>::moveImpl( std::move( other ) );

    other.mNumValuesPerRow = 0;   // force consistency of other with zero storage

    return *this;
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
void ELLStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::StorageWrapper<ELLStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    ContextPtr ctx = getContextPtr();   // will force a valid copy in this context

    if ( other.getFormat() == Format::ELL )
    {
        // both storage have ELL format, use special method for i
        assignELL( static_cast<const ELLStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRData( otherCSR.getNumRows(), otherCSR.getNumColumns(),
                    otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues() );
    }
    else
    {
        HArray<IndexType>  csrIA( ctx );
        HArray<IndexType>  csrJA( ctx );
        HArray<ValueType>  csrValues( ctx );     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        // just a thought for optimization: use mIA, mJA, mValues instead of csrIA, csrJA, csrValues
        // but does not help much at all as resort of entries requires already temporaries.

        setCSRData( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues );
    }
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::assignELL( const ELLStorage<OtherValueType>& other )
{
    if ( static_cast<const _MatrixStorage*>( &other ) == this )
    {
        SCAI_LOG_DEBUG( logger, typeName() << ": self assign, skipped, storage = " << other )
        return;
    }

    SCAI_LOG_DEBUG( logger, "assignELL: other = " << other )

    // The following code is nearly the same as setELLData but takes advantage of given properties
    // i.e. no checks are needed, unused entries in values are not filled up

    auto ctx = getContextPtr();

    // both storage have JDS format, we can just copy the corresponding arrays to the right context

    _MatrixStorage::_assign( other );     // assign member variables of base class

    mNumValuesPerRow = other.getNumValuesPerRow();

    HArrayUtils::assign( mIA, other.getIA(), ctx );
    HArrayUtils::assign( mJA, other.getJA(), ctx );
    HArrayUtils::assign( mValues, other.getValues(), ctx );

    SCAI_LOG_DEBUG( logger, "assignELL: other = " << other << ", this = " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::splitUp(
    IndexType& numRows,
    IndexType& numColumns,
    IndexType& numValuesPerRow,
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values )
{
    numValuesPerRow = mNumValuesPerRow;

    ia = std::move( mIA );
    ja = std::move( mJA );
    values = std::move( mValues );

    mNumValuesPerRow = 0;

    // reset the dimensions of this storage to zero so it remains consistent

    _MatrixStorage::splitUp( numRows, numColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::print( std::ostream& stream ) const
{
    SCAI_LOG_INFO( logger, "print" )
    using std::endl;
    stream << "ELLStorage " << getNumRows() << " x " << getNumColumns() << ", #values = " << getNumValues() << endl;
    ReadAccess<IndexType> ia( mIA );
    ReadAccess<IndexType> ja( mJA );
    ReadAccess<ValueType> values( mValues );

    for ( IndexType i = 0; i < getNumRows(); i++ )
    {
        stream << "Row " << i << " ( " << ia[i] << " entries ) :";

        for ( IndexType jj = 0; jj < ia[i]; ++jj )
        {
            IndexType pos = jj * getNumRows() + i;
            stream << " " << ja[pos] << ":" << values[pos];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format ELLStorage<ValueType>::getFormat() const
{
    return Format::ELL;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType ELLStorage<ValueType>::getNumValues() const
{
    SCAI_LOG_INFO( logger, "getNumValues" )
    IndexType numValues = HArrayUtils::reduce( mIA, common::BinaryOp::ADD, this->getContextPtr() );
    return numValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::purge()
{
    SCAI_LOG_INFO( logger, "purge" )

    _MatrixStorage::setDimension( 0, 0 );

    mNumValuesPerRow = 0;

    mIA.purge();
    mJA.purge();
    mValues.purge();
    mRowIndexes.purge();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_INFO( logger, "set identity # size = " << size )

    _MatrixStorage::setDimension( size, size );

    mNumValuesPerRow = 1;

    const ContextPtr loc = this->getContextPtr();

    // Note: specify template param explicitly as type deduce might fail due to 0, 1

    HArrayUtils::setSameValue<IndexType>( mIA, size, 1, loc );
    HArrayUtils::setSequence<IndexType>( mJA, 0, 1, size, getContextPtr() );
    HArrayUtils::setSameValue<ValueType>( mValues, size, 1, loc );

    SCAI_LOG_INFO( logger, *this << " is identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    _MatrixStorage::setDimension( size, size );

    mNumValuesPerRow = 1;

    // Note: we pass also the current context as it might have been changed

    HArrayUtils::setSameValue<IndexType>( mIA, size, 1, getContextPtr() );
    HArrayUtils::setSequence<IndexType>( mJA, 0, 1, size, getContextPtr() );
    HArrayUtils::assign( mValues, diagonal, getContextPtr() );

    // Note: we do not build row indexes, no row is empty

    SCAI_LOG_INFO( logger, *this << ": diagonal matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::clear()
{
    SCAI_LOG_INFO( logger, "clear" )

    _MatrixStorage::setDimension( 0, 0 );

    mNumValuesPerRow = 0;

    mIA.clear();
    mJA.clear();
    mValues.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::buildCSRSizes( hmemo::HArray<IndexType>& csrSizes ) const
{
    HArrayUtils::assign<IndexType>( csrSizes, mIA );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::buildCSRData(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    _HArray& csrValues ) const
{
    SCAI_REGION( "Storage.ELL.buildCSR" )

    SCAI_LOG_INFO( logger,
                   "buildCSR<" << csrValues.getValueType() << ">"
                   << " from ELL<" << common::getScalarType<ValueType>() << ">"
                   << " on " << *getContextPtr() << " ( preferred )" )

    if ( csrValues.getValueType() == getValueType() )
    {
        auto& castCSRValues = static_cast<HArray<ValueType>&>( csrValues );
        ELLUtils::convertELL2CSR( csrIA, csrJA, castCSRValues, getNumRows(), getNumColumns(), mIA, mJA, mValues, getContextPtr() );
    }
    else
    {
        HArray<ValueType> tmpValues;  // use temporary for conversion of values
        ELLUtils::convertELL2CSR( csrIA, csrJA, tmpValues, getNumRows(), getNumColumns(), mIA, mJA, mValues, getContextPtr() );
        HArrayUtils::_assign( csrValues, tmpValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setCSRData(
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

    if ( csrValues.getValueType() != getValueType() )
    {
        SCAI_LOG_INFO( logger, "setCSRData<" << csrValues.getValueType() << ">, convert values to " << getValueType() )

        HArray<ValueType> sameTypeCSRValues; 
        HArrayUtils::_assign( sameTypeCSRValues, csrValues, getContextPtr() );    
        setCSRData( numRows, numColumns, csrIA, csrJA, sameTypeCSRValues );
        return;
    }
    
    SCAI_LOG_INFO( logger, "setCSRData<" << getValueType() << "> " << numRows << " x " << numColumns << ", nnz = " << csrJA.size() )

    _MatrixStorage::setDimension( numRows, numColumns );

    // csrValues have same type, so we only have to cast

    const auto& sameTypeCSRValues = static_cast<const HArray<ValueType>&>( csrValues );

    ELLUtils::convertCSR2ELL( mIA, mJA, mValues, numRows, numColumns, csrIA, csrJA, sameTypeCSRValues, getContextPtr() );

    if ( numRows > 0 )
    {
        mNumValuesPerRow = mJA.size() / numRows;
    }
    else
    {
        mNumValuesPerRow = 0;
    }

    buildRowIndexes();

    SCAI_LOG_INFO( logger, "ELL: set CSR data done, this = " << *this )

    if ( getNumRows() > 200 && mNumValuesPerRow > 0 )
    {
        // make this check only on larger matrices, dataSize must not be equal 0

        double fillRate = double( csrJA.size() ) / double( mJA.size() );

        if ( fillRate < 0.5 )
        {
            SCAI_LOG_WARN( logger,
                           *this << ": fill rate = " << fillRate << " ( " << csrJA.size() << " non-zero values, "
                           << "but allocated " << mJA.size() << " ), consider using JDS" )
        }
    }
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
    HArrayUtils::assign( mIA, ia, loc );
    HArrayUtils::assign( mJA, ja, loc );

    // _assign must be used here instead of assign as values are untyped
    HArrayUtils::_assign( mValues, values, loc );  // can deal with type conversion

    // fill up my arrays ja and values to make matrix-multiplication fast
    fillValues();

    // check is expensive, so do it only if ASSERT_LEVEL is on DEBUG mode
#ifdef SCAI_ASSERT_LEVEL_DEBUG
    check( "ELLStorage( #row, #cols, #values, #diags, dlg, ilg, perm, ja, values" );
#endif
    SCAI_LOG_INFO( logger, *this << ": set ELLPACK by arrays ia, ja, values" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::fillValues()
{
    static LAMAKernel<ELLKernelTrait::fillELLValues<ValueType> > fillELLValues;
    ContextPtr loc = this->getContextPtr();
    fillELLValues.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ellIA( mIA, loc );
    WriteAccess<IndexType> ellJA( mJA, loc );
    WriteAccess<ValueType> ellValues( mValues, loc );
    fillELLValues[loc]( ellJA.get(), ellValues.get(), ellIA.get(), getNumRows(), mNumValuesPerRow );
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
void ELLStorage<ValueType>::setDiagonal( const ValueType value )
{
    SCAI_LOG_INFO( logger, "setDiagonal for " << *this << ": diagonal value = " << value )

    ELLUtils::setDiagonal( mValues, value, getNumRows(), getNumColumns(), mIA, mJA, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    SCAI_LOG_INFO( logger, "setDiagonalV for " << *this << ": diagonal = " << diagonal )

    SCAI_ASSERT_EQ_ERROR( diagonal.size(), getDiagonalSize(), "diagonal has illegal size" )
    ELLUtils::setDiagonalV( mValues, diagonal, getNumRows(), getNumColumns(), mIA, mJA, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.ELL.getSparseRow" )

    const IndexType nrow  = mIA[i];       // number of non-zero entries in row
    const IndexType offs  = i;            // first non-zero entry
    const IndexType inc   = getNumRows();     // stride between two entries

    // resize the output arrays, invalidate old data before

    jA.clear();
    jA.resize( nrow );
    values.clear();
    values.resize( nrow );

    // just copy the corresponding parts of the csrJA and csrValues array

    common::BinaryOp op = common::BinaryOp::COPY;

    HArrayUtils::setArraySection( jA, 0, 1, mJA, offs, inc, nrow, op, getContextPtr() );
    HArrayUtils::setArraySection( values, 0, 1, mValues, offs, inc, nrow, op, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::getSparseColumn( 
    hmemo::HArray<IndexType>& iA, 
    hmemo::HArray<ValueType>& values, 
    const IndexType j ) const
{   
    SCAI_REGION( "Storage.ELL.getSparseCol" )
    
    // check for legal column index j; but routine works fine and return empty column

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )
    
    HArray<IndexType> columnPositions;  // ellJA[columnPositions[i]] == j for 0 <= i < size

    ELLUtils::getColumnPositions( iA, columnPositions, mIA, mJA, j, getContextPtr() );

    // column values[i] = mValues[ pos[i] ];  

    HArrayUtils::gather( values, mValues, columnPositions, common::BinaryOp::COPY, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::getRow( HArray<ValueType>& row, const IndexType i ) const
{
    SCAI_LOG_TRACE( logger, "getRowImpl # row = " << row << ", i = " << i )
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    static LAMAKernel<ELLKernelTrait::getRow<ValueType> > getRow;
    ContextPtr loc = this->getContextPtr();
    getRow.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    WriteOnlyAccess<ValueType> wRow( row, loc, getNumColumns() );
    const ReadAccess<IndexType> rIa( mIA, loc );
    const ReadAccess<IndexType> rJa( mJA, loc );
    const ReadAccess<ValueType> rValues( mValues, loc );
    getRow[loc]( wRow.get(), i, getNumRows(), getNumColumns(), mNumValuesPerRow, rIa.get(), rJa.get(), rValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.ELL.getCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, getNumRows(), colValues, rowIndexes, ValueType( 0 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const common::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    // ToDo write more efficient kernel routine for setting a row

    ReadAccess<ValueType> rRow( row );

    for ( IndexType j = 0; j < getNumColumns(); ++j )
    {
        if ( rRow[j] == common::Constants::ZERO )
        {
            continue;
        }

        setValue( i, j, static_cast<ValueType>( rRow[j] ), op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const common::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), getNumRows(), "column array to small for set" )

    // ToDo write more efficient kernel routine for setting a column

    ReadAccess<ValueType> rColumn( column );

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        if ( rColumn[i] == common::Constants::ZERO )
        {
            continue;
        }

        setValue( i, j, static_cast<ValueType>( rColumn[i] ), op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    SCAI_LOG_INFO( logger, "getDiagonal for " << *this )
    ELLUtils::getDiagonal( diagonal, getNumRows(), getNumColumns(), mIA, mJA, mValues, getContextPtr() );
    SCAI_ASSERT_EQ_DEBUG( diagonal.size(), getDiagonalSize(), "serious mismatch" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::scale( const ValueType value )
{
    SCAI_LOG_INFO( logger, "scale # value = " << value )
    HArrayUtils::setScalar( mValues, value, common::BinaryOp::MULT, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::UnaryOp::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::scaleRows( const HArray<ValueType>& values )
{
    // (MULT)iply each row with an individual value

    ELLUtils::setRows( mValues, getNumRows(), getNumColumns(), mIA, mJA, values, common::BinaryOp::MULT, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::scaleColumns( const HArray<ValueType>& values )
{
    // (MULT)iply each column with an individual value

    ELLUtils::setColumns( mValues, getNumRows(), getNumColumns(), mIA, mJA, values, common::BinaryOp::MULT, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& ELLStorage<ValueType>::getIA() const
{
    return mIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& ELLStorage<ValueType>::getJA() const
{
    return mJA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& ELLStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::~ELLStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~ELLStorage for matrix " << getNumRows() << " x " << getNumColumns() << ", # nnr = " << mNumValuesPerRow )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::check( const char* msg ) const
{
    SCAI_LOG_INFO( logger, "check # msg = " << msg )
    SCAI_ASSERT_EQ_ERROR( getNumRows(), mIA.size(), "inconsistent ELLStorage: " << *this )
    SCAI_ASSERT_EQ_ERROR( mNumValuesPerRow * getNumRows(), mJA.size(), "inconsistent ELLStorage: " << *this )
    SCAI_ASSERT_EQ_ERROR( mJA.size(), mValues.size(), "inconsistent ELLStorage: " << *this )
    static LAMAKernel<ELLKernelTrait::check> check;
    ContextPtr loc = this->getContextPtr();
    check.getSupportedContext( loc );
    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );
    SCAI_CONTEXT_ACCESS( loc )
    check[loc]( getNumRows(), mNumValuesPerRow, getNumColumns(), rIa.get(), rJa.get(), msg );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate ELL sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    _MatrixStorage::setDimension( numRows, numColumns );

    mIA.clear();
    mIA.resize( getNumRows() );
    HArrayUtils::setScalar( mIA, IndexType( 0 ), common::BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "ELLStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << getNumRows() << " x " << getNumColumns()
           << ", nnr = " << mNumValuesPerRow << ", threshold = " << mCompressThreshold 
           << ", ia = " << mIA << ", ja = " << mJA << ", values = " << mValues << " )"; 
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    IndexType pos = ELLUtils::getValuePos( i, j, mIA, mJA, getContextPtr() );

    ValueType val = 0;

    if ( pos != invalidIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, getNumRows() * mNumValuesPerRow,
                                       "illegal value position for ( " << i << ", " << j << " )" );

        val = mValues[ pos ];
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const common::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_DEBUG( logger, "set value (" << i << ", " << j << ")" )

    IndexType pos = ELLUtils::getValuePos( i, j, mIA, mJA, getContextPtr() );

    if ( pos == invalidIndex )
    {
        COMMON_THROWEXCEPTION( "ELL storage has no entry ( " << i << ", " << j << " ) " )
    }

    utilskernel::HArrayUtils::setVal( mValues, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::prefetch( const ContextPtr context ) const
{
    SCAI_LOG_INFO( logger, "prefetch # context " << context )
    SCAI_LOG_DEBUG( logger, "Starting prefetch of " << *this << " to " << context )
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
void ELLStorage<ValueType>::buildRowIndexes()
{
    // build row indexes if there are only few rows that are not empty (e.g. for Halo storage)

    ELLUtils::nonEmptyRows( mRowIndexes, mIA, mCompressThreshold, getContextPtr() );

    SCAI_LOG_INFO( logger, "#row indexes = " << mRowIndexes.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::compress( const RealType<ValueType> eps )
{
    SCAI_LOG_INFO( logger, "compress: eps = " << eps )

    ELLUtils::compress( mIA, mJA, mValues, mNumValuesPerRow, eps, getContextPtr() );

    buildRowIndexes();   // sizes of rows might have changed
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::swap( ELLStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "swap # other = " << other )

    // swap base class

    MatrixStorage<ValueType>::swap( other );

    // swap my member variables

    std::swap( mNumValuesPerRow, other.mNumValuesPerRow );
    mIA.swap( other.mIA );
    mJA.swap( other.mJA );
    mValues.swap( other.mValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t ELLStorage<ValueType>::getMemoryUsageImpl() const
{
    SCAI_LOG_INFO( logger, "getMemoryUsageImpl" )
    size_t memoryUsage = 0;
    memoryUsage += sizeof( IndexType );
    memoryUsage += sizeof( IndexType ) * mIA.size();
    memoryUsage += sizeof( IndexType ) * mJA.size();
    memoryUsage += sizeof( ValueType ) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::matrixTimesVector(
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

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::matrixTimesVectorAsync(
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

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_REGION( "Storage.ELL.gemv" )

    SCAI_LOG_INFO( logger,
                   "gemv<" << getValueType() << "> ( op = " << op << ", async = " << async 
                   << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    MatrixStorage<ValueType>::gemvCheck( alpha, x, beta, y, op );  // checks for correct sizes

    SyncToken* token = NULL;

    if ( beta == common::Constants::ZERO )
    {
        // take version that does not access y at all (can be undefined or aliased to result)

        token = ELLUtils::gemv0( result, alpha, x,
                                 getNumRows(), getNumColumns(), mNumValuesPerRow, mIA, mJA, mValues,
                                 op, async, getContextPtr() );
    }
    else if ( &result == &y && ( beta == common::Constants::ONE ) && ( mRowIndexes.size() > 0 ) )
    {
        // y += A * x,  where only some rows in A are filled, uses more efficient routine

        token = ELLUtils::gemvSp( result, alpha, x, getNumRows(), getNumColumns(), mNumValuesPerRow,
                                  mIA, mJA, mValues, op, mRowIndexes, async, getContextPtr() );
    }
    else
    {
        token = ELLUtils::gemv( result, alpha, x, beta, y,
                                getNumRows(), getNumColumns(), mNumValuesPerRow, 
                                mIA, mJA, mValues,
                                op, async, getContextPtr() );
    }

    return token;
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

    bool async = false;  // no sync token, call will return NULL

    ELLUtils::jacobi( solution, omega, oldSolution, rhs, mIA, mJA, mValues, async, getContextPtr() );
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

    bool async = true;  // call will return valid SyncToken

    SyncToken* token = ELLUtils::jacobi( solution, omega, oldSolution, rhs, mIA, mJA, mValues, async, getContextPtr() );

    if ( token == NULL )
    {
        // there was no asynchronous execution at all

        token = new NoSyncToken();
    }
 
    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.ELL.jacobiIterateHalo" )

    SCAI_LOG_INFO( logger, "HOST: Jacobi iteration on halo matrix data." )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumColumns(), oldHaloSolution.size() )

    ELLUtils::jacobiHalo( localSolution, omega, localDiagonal, oldHaloSolution,
                          mIA, mJA, mValues, mRowIndexes, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{   
    halo.halo2Global( mJA );
    _MatrixStorage::setDimension( getNumRows(), globalNumColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> ELLStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    return HArrayUtils::l1Norm( mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> ELLStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    // Note: un-used entries of values have been filled with 0, so use norm for arrays

    return HArrayUtils::l2Norm( mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> ELLStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    return HArrayUtils::maxNorm( mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const ValueType beta,
    const MatrixStorage<ValueType>& b )
{
    SCAI_LOG_INFO( logger, "this = " << alpha << " * A + " << beta << " * B" << ", with A = " << a << ", B = " << b )

    SCAI_REGION( "Storage.ELL.plusMatrix" )

    if ( a.getFormat() == Format::ELL )
    {
        const ELLStorage<ValueType>& ellA = static_cast<const ELLStorage<ValueType>&>( a );

        if ( b.getFormat() == Format::ELL )
        {   
            matrixAddMatrixELL( alpha, ellA, beta, static_cast<const ELLStorage<ValueType>&>( b ) );
        }
        else
        {   
            matrixAddMatrixELL( alpha, ellA, beta, convert<ELLStorage<ValueType>>( b ) );
        }
    }
    else
    {
        auto ellA = convert<ELLStorage<ValueType>>( a );
   
        if ( b.getFormat() == Format::ELL )
        {
            matrixAddMatrixELL( alpha, ellA, beta, static_cast<const ELLStorage<ValueType>&>( b ) );
        }
        else
        {
            matrixAddMatrixELL( alpha, ellA, beta, convert<ELLStorage<ValueType>>( b ) );
        }
    }
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
    SCAI_REGION( "Storage.ELL.timesMatrix" )

    SCAI_LOG_INFO( logger,
                   "this = " << alpha << " * A * B + " << beta << " * C, with " << "A = " << a << ", B = " << b << ", C = " << c )
    const ELLStorage<ValueType>* ellA = NULL;
    const ELLStorage<ValueType>* ellB = NULL;
    const ELLStorage<ValueType>* ellC = NULL;
    //    std::shared_ptr<ELLStorage<ValueType> > tmpA;
    //    std::shared_ptr<ELLStorage<ValueType> > tmpB;
    std::shared_ptr<ELLStorage<ValueType> > tmpC;

    if ( a.getFormat() == Format::ELL )
    {
        ellA = dynamic_cast<const ELLStorage<ValueType>*>( &a );
        SCAI_ASSERT_DEBUG( ellA, "could not cast to ELLStorage " << a )
    }
    else
    {
        SCAI_LOG_ERROR( logger, a << ": a not ELL format" )
    }

    if ( b.getFormat() == Format::ELL )
    {
        ellB = dynamic_cast<const ELLStorage<ValueType>*>( &b );
        SCAI_ASSERT_DEBUG( ellB, "could not cast to ELLStorage " << b )
    }
    else
    {
        SCAI_UNSUPPORTED( b << ": b not ELL format" )
    }

    if ( ellA == NULL || ellB == NULL )
    {
        // input matrices not ELL format, so try via ELL
        MatrixStorage<ValueType>::matrixTimesMatrix( alpha, a, b, beta, c );
        return;
    }

    if ( beta != scai::common::Constants::ZERO )
    {
        if ( ( c.getFormat() == Format::ELL ) && ( &c != this ) )
        {
            ellC = dynamic_cast<const ELLStorage<ValueType>*>( &c );
            SCAI_ASSERT_DEBUG( ellC, "could not cast to ELLStorage " << c )
        }
        else
        {
            SCAI_UNSUPPORTED( c << ": ELL temporary required for matrix add" )
            tmpC.reset( new ELLStorage<ValueType>() );
            tmpC->assign( c );
            ellC = tmpC.get();
        }
    }

    ELLStorage<ValueType> tmp;
    tmp.matrixTimesMatrixELL( alpha, *ellA, *ellB );

    if ( beta != scai::common::Constants::ZERO )
    {
        ELLStorage<ValueType> tmp1;
        tmp1.matrixAddMatrixELL( static_cast<ValueType>( 1.0 ), tmp, beta, *ellC );
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

    if ( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary

        ELLStorage<ValueType> tmp;
        tmp.matrixTimesMatrixELL( alpha, a, b );
        swap( tmp ); // safe as tmp will be destroyed afterwards
        return;
    }

    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<ELLKernelTrait::matrixMultiplySizes> matrixMultiplySizes;
    static LAMAKernel<ELLKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    ContextPtr loc = Context::getHostPtr();  // not yet available on other devices

    SCAI_ASSERT_ERROR( &a != this, "matrixTimesMatrix: alias of a with this result matrix" )
    SCAI_ASSERT_ERROR( &b != this, "matrixTimesMatrix: alias of b with this result matrix" )
    SCAI_ASSERT_EQUAL_ERROR( a.getNumColumns(), b.getNumRows() )
    allocate( a.getNumRows(), b.getNumColumns() );
    {
        ReadAccess<IndexType> aIA( a.getIA(), loc );
        ReadAccess<IndexType> aJA( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );
        ReadAccess<IndexType> bIA( b.getIA(), loc );
        ReadAccess<IndexType> bJA( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );
        allocate( a.getNumRows(), b.getNumColumns() );
        WriteOnlyAccess<IndexType> cIA( mIA, loc, getNumRows() );
        SCAI_CONTEXT_ACCESS( loc )
        // 1. Step: compute resulting IA array
        matrixMultiplySizes[loc] ( cIA.get(), a.getNumRows(), a.getNumColumns(), b.getNumRows(), false, aIA.get(), aJA.get(),
                                   a.getNumValuesPerRow(), bIA.get(), bJA.get(), b.getNumValuesPerRow() );
        // 2. Step: compute length of longest row
        mNumValuesPerRow = reduce[ loc ]( cIA.get(), getNumRows(), 0, common::BinaryOp::MAX );
        // 3. Step: Allocate IA and Values arrays with new size
        WriteOnlyAccess<IndexType> cJA( mJA, loc, mNumValuesPerRow * getNumRows() );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValuesPerRow * getNumRows() );
        // 4. Step: Compute cJA and cValues
        matrixMultiply[loc]( cJA.get(), cValues.get(), cIA.get(), mNumValuesPerRow, getNumRows(), getNumColumns(), b.getNumRows(),
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

    if ( &a == this || &b == this )
    {
        // due to alias we would get problems with Write/Read access, so use a temporary
        ELLStorage<ValueType> tmp;
        tmp.setContextPtr( this->getContextPtr() );
        tmp.matrixAddMatrixELL( alpha, a, beta, b );
        swap( tmp ); // safe as tmp will be destroyed afterwards
        return;
    }

    static LAMAKernel<ELLKernelTrait::matrixAddSizes> matrixAddSizes;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    static LAMAKernel<ELLKernelTrait::matrixAdd<ValueType> > matrixAdd;
    ContextPtr loc = this->getContextPtr();
    matrixAddSizes.getSupportedContext( loc, reduce, matrixAdd );
    SCAI_ASSERT_ERROR( &a != this, "matrixAddMatrix: alias of a with this result matrix" )
    SCAI_ASSERT_ERROR( &b != this, "matrixAddMatrix: alias of b with this result matrix" )
    allocate( a.getNumRows(), a.getNumColumns() );
    SCAI_ASSERT_EQUAL_ERROR( getNumRows(), b.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), b.getNumColumns() )
    {
        ReadAccess<IndexType> aIA( a.getIA(), loc );
        ReadAccess<IndexType> aJA( a.getJA(), loc );
        ReadAccess<ValueType> aValues( a.getValues(), loc );
        ReadAccess<IndexType> bIA( b.getIA(), loc );
        ReadAccess<IndexType> bJA( b.getJA(), loc );
        ReadAccess<ValueType> bValues( b.getValues(), loc );
        WriteOnlyAccess<IndexType> cIA( mIA, loc, getNumRows() );
        SCAI_CONTEXT_ACCESS( loc )
        // 1. Step: Compute IA array
        matrixAddSizes[loc]( cIA.get(), a.getNumRows(), a.getNumColumns(), false, aIA.get(), aJA.get(),
                             a.getNumValuesPerRow(), bIA.get(), bJA.get(), b.getNumValuesPerRow() );
        // 2. Step: compute length of longest row
        mNumValuesPerRow = reduce[loc]( cIA.get(), getNumRows(), 0, common::BinaryOp::MAX );
        // 3. Step: Allocate IA and Values arrays with new size
        WriteOnlyAccess<IndexType> cJA( mJA, loc, mNumValuesPerRow * getNumRows() );
        WriteOnlyAccess<ValueType> cValues( mValues, loc, mNumValuesPerRow * getNumRows() );
        // 4. Step: Compute cJA and cValues
        matrixAdd[loc]( cJA.get(), cValues.get(), cIA.get(), mNumValuesPerRow, getNumRows(), getNumColumns(), false, alpha,
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
ELLStorage<ValueType>* ELLStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    std::unique_ptr<ELLStorage<ValueType> > storage( new ELLStorage<ValueType>( getContextPtr() ) );
    storage->allocate( numRows, numColumns );
    return storage.release();
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

template<typename ValueType>
std::string ELLStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "ELLStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* ELLStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
const char* ELLStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType ELLStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::ELL, common::getScalarType<ValueType>() );
}

template<typename ValueType>
MatrixStorageCreateKeyType ELLStorage<ValueType>::getCreateValue() const
{
    return createValue();
}

template<typename ValueType>
_MatrixStorage* ELLStorage<ValueType>::create()
{
    return new ELLStorage<ValueType>();
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( ELLStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
