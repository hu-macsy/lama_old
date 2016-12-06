/**
 * @file ELLStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/instantiate.hpp>

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
    ContextPtr context )

    : CRTPMatrixStorage<ELLStorage<ValueType>, ValueType>( numRows, numColumns ), mNumValuesPerRow( 0 )
{
    setContextPtr( context );
    // Initialization requires correct values for the IA array with 0
    mIA.resize( mNumRows );
    // ellSizes[] = 0 @ context
    HArrayUtils::setScalar( mIA, IndexType( 0 ), utilskernel::binary::COPY, context );
    SCAI_LOG_DEBUG( logger, "ELLStorage for matrix " << mNumRows << " x " << mNumColumns << ", no elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage()

    : CRTPMatrixStorage<ELLStorage<ValueType>, ValueType>( 0, 0 ), mNumValuesPerRow( 0 )
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

    : CRTPMatrixStorage<ELLStorage<ValueType>, ValueType>()
{
    SCAI_LOG_INFO( logger, "constructor with ELL data array" )
    setELLData( numRows, numColumns, numValuesPerRows, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ELLStorage<ValueType>::ELLStorage( const ELLStorage<ValueType>& other )

    : CRTPMatrixStorage<ELLStorage<ValueType>, ValueType>( 0, 0 )
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
    if ( &other == this )
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
void ELLStorage<ValueType>::print( std::ostream& stream ) const
{
    SCAI_LOG_INFO( logger, "print" )
    using std::endl;
    stream << "ELLStorage " << mNumRows << " x " << mNumColumns << ", #values = " << getNumValues() << endl;
    ReadAccess<IndexType> ia( mIA );
    ReadAccess<IndexType> ja( mJA );
    ReadAccess<ValueType> values( mValues );

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        stream << "Row " << i << " ( " << ia[i] << " entries ) :";

        for ( IndexType jj = 0; jj < ia[i]; ++jj )
        {
            IndexType pos = jj * mNumRows + i;
            stream << " " << ja[pos] << ":" << values[pos];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format::MatrixStorageFormat ELLStorage<ValueType>::getFormat() const
{
    return Format::ELL;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType ELLStorage<ValueType>::getNumValues() const
{
    SCAI_LOG_INFO( logger, "getNumValues" )
    IndexType numValues = HArrayUtils::reduce( mIA, utilskernel::binary::ADD, this->getContextPtr() );
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
    const ContextPtr loc = this->getContextPtr();
    mIA.clear();
    mIA.resize( mNumRows );
    HArrayUtils::setScalar( mIA, IndexType( 1 ), utilskernel::binary::COPY, loc );
    HArrayUtils::setOrder( mJA, mNumRows );
    mValues.clear();
    mValues.resize( mNumRows );
    HArrayUtils::setScalar( mValues, ValueType( 1 ), utilskernel::binary::COPY, loc );
    mDiagonalProperty = true;
    SCAI_LOG_INFO( logger, *this << " is identity matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool ELLStorage<ValueType>::checkDiagonalProperty() const
{
    SCAI_LOG_INFO( logger, "checkDiagonalProperty" )
    IndexType numDiagonals = common::Math::min( mNumRows, mNumColumns );
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
        static LAMAKernel<ELLKernelTrait::hasDiagonalProperty> ellHasDiagonalProperty;
        // check it where the JA array has a valid copy
        ContextPtr loc = mJA.getValidContext();
        ellHasDiagonalProperty.getSupportedContext( loc );
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
    SCAI_REGION( "Storage.ELL.buildCSR" )
    SCAI_LOG_INFO( logger,
                   "buildCSR<" << common::getScalarType<OtherValueType>() << ">"
                   << " from ELL<" << common::getScalarType<ValueType>() << ">"
                   << " on " << *context << " ( preferred )" )
    // step 1 : compute IA offsets
    IndexType numValues = 0;
    ia.clear();
    ia.reserve( context, mNumRows + 1 );  // reserve one more entry
    HArrayUtils::setArrayImpl( ia, mIA, utilskernel::binary::COPY, context );

    if ( ja == NULL || values == NULL )
    {
        return;
    }

    numValues = HArrayUtils::scan( ia, context );
    // step 2 : compute the arrays ja and values
    static LAMAKernel<ELLKernelTrait::getCSRValues<ValueType, OtherValueType> > getCSRValues;
    ContextPtr loc = context;
    getCSRValues.getSupportedContext( loc );
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
    SCAI_REGION( "Storage.ELL.setCSR" )
    SCAI_LOG_INFO( logger,
                   "set CSR data on " << *context << ": numRows = " << numRows << ", numColumns = " << numColumns
                   << ", numValues = " << numValues << ", compress threshold = " << mCompressThreshold )

    if ( numRows == 0 )
    {
        // just allocate will clear member arrays
        allocate( numRows, numColumns );
        return;
    }

    _MatrixStorage::setDimension( numRows, numColumns );
    // build array with non-zero values per row
    common::unique_ptr<HArray<IndexType> > tmpOffsets;
    const HArray<IndexType>* offsets = &ia;

    if ( ia.size() == numRows + 1 )
    {
        ContextPtr loc = context;
        static LAMAKernel<CSRKernelTrait::offsets2sizes > offsets2sizes;
        offsets2sizes.getSupportedContext( loc );
        ReadAccess<IndexType> csrIA( ia, loc );
        WriteOnlyAccess<IndexType> ellSizes( mIA, loc, mNumRows );
        SCAI_CONTEXT_ACCESS( loc )
        offsets2sizes[ loc ]( ellSizes.get(), csrIA.get(), mNumRows );
    }
    else if ( ia.size() == numRows )
    {
        HArrayUtils::assign( mIA, ia, context );
        // as the offset array is also needed
        tmpOffsets.reset( ia.copy() );
        IndexType total = HArrayUtils::scan( *tmpOffsets, context );
        SCAI_ASSERT_EQUAL( total, numValues, "sizes do not sum up correctly" )
        offsets = tmpOffsets.get();
    }
    else
    {
        COMMON_THROWEXCEPTION( "ia array has illegal size " << ia.size() << " for #rows = " << numRows )
    }

    // determine the maximal number of non-zero in one row
    mNumValuesPerRow = HArrayUtils::reduce( mIA, utilskernel::binary::MAX, context );
    SCAI_LOG_DEBUG( logger, "setCSRData, #values/row = " << mNumValuesPerRow )
    //  Now we know the size of the ja and values arrays for the ELL format
    const IndexType dataSize = mNumValuesPerRow * mNumRows;

    if ( mNumRows > 200 && mNumValuesPerRow > 0 )
    {
        // make this check only on larger matrices, dataSize must not be equal 0
        double fillRate = double( numValues ) / double( dataSize );

        if ( fillRate < 0.5 )
        {
            SCAI_LOG_WARN( logger,
                           *this << ": fill rate = " << fillRate << " ( " << numValues << " non-zero values ), consider using JDS" )
        }
    }

    // Get function pointers for needed routines at the LAMA interface
    static LAMAKernel<ELLKernelTrait::hasDiagonalProperty > hasDiagonalProperty;
    static LAMAKernel<ELLKernelTrait::setCSRValues<ValueType, OtherValueType> > setCSRValues;
    ContextPtr loc = context;
    setCSRValues.getSupportedContext( loc, hasDiagonalProperty );
    {
        // now fill the matrix values and column indexes
        ReadAccess<IndexType> csrIA( *offsets, loc );
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
        SCAI_LOG_DEBUG( logger, " size = " << ellJA.size() )
        IndexType numDiagonals = std::min( mNumRows, mNumColumns );

        if ( numDiagonals == 0 )
        {
            mDiagonalProperty = true;
        }
        else if ( numValues == 0 )
        {
            mDiagonalProperty = false;
        }
        else
        {
            SCAI_CONTEXT_ACCESS( loc )
            mDiagonalProperty = hasDiagonalProperty[loc]( numDiagonals, ellJA.get() );
        }
    }

    if ( numRows == numColumns && !mDiagonalProperty )
    {
        SCAI_LOG_INFO( logger, *this << ": square matrix has not diagonal property" )
    }

    buildRowIndexes( loc );
    SCAI_LOG_DEBUG( logger, "convert CSR -> ELL done: " << *this )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::setDIADataImpl(
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
    HArrayUtils::setArrayImpl( mIA, ia, utilskernel::binary::COPY, loc );
    HArrayUtils::setArrayImpl( mJA, ja, utilskernel::binary::COPY, loc );
    // setArray must be used here instead of setArrayImpl as values is untyped 
    HArrayUtils::setArray( mValues, values, utilskernel::binary::COPY, loc );  // also type conversion
    // fill up my arrays ja and values to make matrix-multiplication fast
    {
        static LAMAKernel<ELLKernelTrait::fillELLValues<ValueType> > fillELLValues;
        ContextPtr loc = this->getContextPtr();
        fillELLValues.getSupportedContext( loc );
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
    ContextPtr loc = this->getContextPtr();
    setVal.getSupportedContext( loc );
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( mValues, loc );
    setVal[ loc ]( wValues.get(), numDiagonalElements, value, utilskernel::binary::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::setDiagonalImpl( const HArray<OtherType>& diagonal )
{
    SCAI_LOG_INFO( logger, "setDiagonalImpl # diagonal = " << diagonal )
    IndexType numDiagonalElements = std::min( mNumColumns, mNumRows );
    static LAMAKernel<UtilKernelTrait::set<ValueType, OtherType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    // ELL format with diagonal property: diagonal is just the first column in mValues
    set[ loc ]( wValues.get(), rDiagonal.get(), numDiagonalElements, utilskernel::binary::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_LOG_TRACE( logger, "getRowImpl # row = " << row << ", i = " << i )
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    static LAMAKernel<ELLKernelTrait::getRow<ValueType, OtherType> > getRow;
    ContextPtr loc = this->getContextPtr();
    getRow.getSupportedContext( loc );
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
void ELLStorage<ValueType>::getColumnImpl( HArray<OtherType>& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.ELL.getCol" )

    static LAMAKernel<ELLKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<IndexType> valuePos;     // positions in the values array
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, mNumRows );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumRows );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), mNumRows, rJA.get(), mNumValuesPerRow );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    column.init( ValueType( 0 ), mNumRows );

    // column[ row ] = mValues[ pos ];

    HArrayUtils::gatherImpl( colValues, mValues, valuePos, utilskernel::binary::COPY, loc );
    HArrayUtils::scatterImpl( column, rowIndexes, colValues, utilskernel::binary::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::setRowImpl( const HArray<OtherType>& row, const IndexType i,
                                        const utilskernel::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), mNumColumns, "row array to small for set" )

    // ToDo write more efficient kernel routine for setting a row

    ReadAccess<OtherType> rRow( row );

    for ( IndexType j = 0; j < mNumColumns; ++j )
    {
        if ( rRow[j] == common::constants::ZERO )
        {
            continue;
        }

        setValue( i, j, static_cast<ValueType>( rRow[j] ), op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::setColumnImpl( const HArray<OtherType>& column, const IndexType j,
                                           const utilskernel::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), mNumRows, "column array to small for set" )

    // ToDo write more efficient kernel routine for setting a column

    ReadAccess<OtherType> rColumn( column );

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        if ( rColumn[i] == common::constants::ZERO )
        {
            continue;
        }

        setValue( i, j, static_cast<ValueType>( rColumn[i] ), op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void ELLStorage<ValueType>::getDiagonalImpl( HArray<OtherType>& diagonal ) const
{
    SCAI_LOG_INFO( logger, "getDiagonalImpl # diagonal = " << diagonal )
    IndexType numDiagonalElements = common::Math::min( mNumColumns, mNumRows );
    // OtherType is output type, so use it as first template argument
    static LAMAKernel<UtilKernelTrait::set<OtherType, ValueType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );
    // ELL format with diagonal property: diagonal is just the first column in mValues
    SCAI_CONTEXT_ACCESS( loc )
    set[loc]( wDiagonal.get(), rValues.get(), numDiagonalElements, utilskernel::binary::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::scaleImpl( const ValueType value )
{
    SCAI_LOG_INFO( logger, "scaleImpl # value = " << value )
    HArrayUtils::setScalar( mValues, value, utilskernel::binary::MULT, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, utilskernel::unary::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::scaleImpl( const HArray<OtherValueType>& values )
{
    SCAI_LOG_INFO( logger, "scaleImpl # values = " << values )
    static LAMAKernel<ELLKernelTrait::scaleValue<ValueType, OtherValueType> > scaleValue;
    ContextPtr loc = this->getContextPtr();
    scaleValue.getSupportedContext( loc );
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
    ContextPtr loc = this->getContextPtr();
    check.getSupportedContext( loc );
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
        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );
        SCAI_CONTEXT_ACCESS( loc )
        WriteOnlyAccess<IndexType> ia( mIA, loc, mNumRows );
        setVal[ loc ]( ia.get(), mNumRows, 0, utilskernel::binary::COPY );
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
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<ELLKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, mNumRows, mNumValuesPerRow, rIa.get(), rJa.get() );

    ValueType val = 0;

    if ( pos != nIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumRows * mNumValuesPerRow,
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
                                      const utilskernel::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_DEBUG( logger, "set value (" << i << ", " << j << ")" )

    static LAMAKernel<ELLKernelTrait::getValuePos> getValuePos;
    
    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    
    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, mNumRows, mNumValuesPerRow, rIa.get(), rJa.get() );

    if ( pos == nIndex )
    {
        COMMON_THROWEXCEPTION( "ELL storage has no entry ( " << i << ", " << j << " ) " )
    }

    utilskernel::HArrayUtils::setValImpl( mValues, pos, val, op );
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
void ELLStorage<ValueType>::buildRowIndexes( const ContextPtr context )
{
    SCAI_LOG_INFO( logger, "buildRowIndexes # loc = " << context )
    mRowIndexes.clear();

    if ( mNumRows == 0 )
    {
        return;
    }

    // Note: compress functionality in HArrayUtils available but we
    // reimplement it here in the same way as compress is optionally done
    // depending on the threshold value

    // Get function pointers for needed kernel routines

    static LAMAKernel<UtilKernelTrait::countNonZeros<IndexType> > countNonZeros;
    static LAMAKernel<UtilKernelTrait::compress<IndexType> > compress;

    // choose location where both routines are available

    ContextPtr loc = context;
    countNonZeros.getSupportedContext( loc, compress );

    ReadAccess<IndexType> ellIA( mIA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // count the number of non-zero rows to have a good value for allocation of rowIndexes

    IndexType nonZeroRows = countNonZeros[loc]( ellIA.get(), mNumRows, 0 );

    float usage = float( nonZeroRows ) / float( mNumRows );

    if ( usage >= mCompressThreshold )
    {
        SCAI_LOG_INFO( logger,
                       "ELLStorage: do not build row indexes, usage = " << usage << " >= " << mCompressThreshold << " ( threshold )" )
        return;
    }

    WriteOnlyAccess<IndexType> rowIndexes( mRowIndexes, loc, nonZeroRows );

    IndexType cnt = compress[loc]( NULL, rowIndexes.get(), ellIA.get(), mNumRows, 0 );

    SCAI_ASSERT_EQ_ERROR( cnt, nonZeroRows, "serious mismatch" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::compress( const ValueType eps /* = 0.0 */ )
{
    SCAI_LOG_INFO( logger, "compress: eps = " << eps )

    ContextPtr loc = this->getContextPtr();
    static LAMAKernel<ELLKernelTrait::compressIA<ValueType> > compressIA;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    compressIA.getSupportedContext( loc, reduce );

    IndexType newNumValuesPerRow = -1;

    LArray<IndexType> newIAArray;
    {
        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> IA( mIA, loc );
        ReadAccess<IndexType> JA( mJA, loc );
        ReadAccess<ValueType> values( mValues, loc );
        // 1. Step: Check for 0 elements and write new IA array
        WriteOnlyAccess<IndexType> newIA( newIAArray, loc, mNumRows );
        compressIA[loc]( IA.get(), JA.get(), values.get(), mNumRows, mNumValuesPerRow, eps, newIA.get() );
        // 2. Step: compute length of longest row
        newNumValuesPerRow = reduce[ loc ]( newIA.get(), mNumRows, 0, utilskernel::binary::MAX );
    }

    // Do further steps, if new array could be smaller
    if ( newNumValuesPerRow < mNumValuesPerRow )
    {
        static LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;
        compressValues.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        // 3. Step: Allocate new JA and Values array
        LArray<ValueType> newValuesArray;
        LArray<IndexType> newJAArray;

        {
            ReadAccess<IndexType> IA( mIA, loc );
            ReadAccess<IndexType> JA( mJA, loc );
            ReadAccess<ValueType> values( mValues, loc );
            WriteOnlyAccess<ValueType> newValues( newValuesArray, loc, mNumRows * newNumValuesPerRow );
            WriteOnlyAccess<IndexType> newJA( newJAArray, loc, mNumRows * newNumValuesPerRow );
            // 4. Step: Compute new JA and Values array
            compressValues[loc]( IA.get(), JA.get(), values.get(), mNumRows, mNumValuesPerRow, eps, newNumValuesPerRow,
                                 newJA.get(), newValues.get() );
        }

        mIA.swap( newIAArray );
        mJA.swap( newJAArray );
        mValues.swap( newValuesArray );
        mNumValuesPerRow = newNumValuesPerRow;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::swap( _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( getFormat(), other.getFormat(), "swap only for same storage format" )
    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(), "swap only for same value type" )

    // only in debug mode use the more expensive dynamic cast for verification

    SCAI_ASSERT_DEBUG( dynamic_cast<ELLStorage<ValueType>* >( &other ), "illegal storage to swap" )

    swapImpl( reinterpret_cast<ELLStorage<ValueType>& >( other ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::swapImpl( ELLStorage<ValueType>& other )
{
    SCAI_LOG_INFO( logger, "swap # other = " << other )

    MatrixStorage<ValueType>::swapMS( other );

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
    bool async = false; // synchronously execution, no SyncToken required
    SyncToken* token = gevm( result, alpha, x, beta, y, async );
    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
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
        HArrayUtils::binaryOpScalar1( result, beta, y, utilskernel::binary::MULT, this->getContextPtr() );

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
SyncToken* ELLStorage<ValueType>::gevm(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool  async ) const
{
    SCAI_REGION( "Storage.ELL.gevm" )
    SCAI_LOG_INFO( logger,
                   "GEVM ( async = " << async << " ), result = " << alpha << " * A * x + " << beta << " * y "
                   << ", result = " << result << ", x = " << x << ", y = " << y
                   << ", A (this) = " << *this );

    if ( alpha == common::constants::ZERO || ( mNumValuesPerRow == 0 ) )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::binaryOpScalar1( result, beta, y, utilskernel::binary::MULT, this->getContextPtr() );

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
        // y += A * x,  where only some rows in A are filled, uses more efficient routine
        return sparseGEVM( result, alpha, x, async );
    }
    else
    {
        return normalGEVM( result, alpha, x, beta, y, async );
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
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_CONTEXT_ACCESS( loc )
    SCAI_ASYNCHRONOUS( syncToken.get() )
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
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::normalGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::normalGEVM<ValueType> > normalGEVM;
    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() )
    SCAI_CONTEXT_ACCESS( loc )
    // Note: alias &result == &y possible
    //       ReadAccess on y before WriteOnlyAccess on result guarantees valid data
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
    normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(),
                     mNumRows, mNumColumns, mNumValuesPerRow,
                     ellIA.get(), ellJA.get(), ellValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( ellIA.releaseDelayed() );
        syncToken->pushRoutine( ellJA.releaseDelayed() );
        syncToken->pushRoutine( ellValues.releaseDelayed() );
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
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() )
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
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::normalGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::normalGEVM<ValueType> > normalGEVM;
    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() )
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, mNumColumns );
    normalGEVM[loc]( wResult.get(), alpha, rX.get(), 0, NULL, mNumRows, mNumColumns, mNumValuesPerRow,
                     ellIA.get(), ellJA.get(), ellValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( ellIA.releaseDelayed() );
        syncToken->pushRoutine( ellJA.releaseDelayed() );
        syncToken->pushRoutine( ellValues.releaseDelayed() );
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
    ContextPtr loc = this->getContextPtr();
    sparseGEMV.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() )
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
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLStorage<ValueType>::sparseGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::sparseGEVM<ValueType> > sparseGEVM;
    ContextPtr loc = this->getContextPtr();
    sparseGEVM.getSupportedContext( loc );
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() )
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );
    // result += alpha * thisMatrix * x, can take advantage of row indexes
    IndexType numNonZeroRows = mRowIndexes.size();
    ReadAccess<IndexType> rRowIndexes( mRowIndexes, loc );
    sparseGEVM[loc]( wResult.get(), alpha, rX.get(),
                     mNumRows, mNumColumns, mNumValuesPerRow,
                     numNonZeroRows, rRowIndexes.get(),
                     ellIA.get(), ellJA.get(), ellValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( rRowIndexes.releaseDelayed() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( ellIA.releaseDelayed() );
        syncToken->pushRoutine( ellJA.releaseDelayed() );
        syncToken->pushRoutine( ellValues.releaseDelayed() );
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
    bool async = true;
    SyncToken* token = gevm( result, alpha, x, beta, y, async );
    SCAI_ASSERT( token, "NULL token not allowed for asynchronous execution gemv, alpha = " << alpha << ", beta = " << beta )
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
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": jacobiIterate requires diagonal property" )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    // matrix must be square
    static LAMAKernel<ELLKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    // make all needed data available at loc
    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );
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
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );

    if ( loc->getType() == Context::Host )
    {
        // used later in OpenMP to generate a TaskSyncToken
        void ( ELLStorage::*jb )(
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
    SCAI_ASSERT_EQUAL_DEBUG( mNumRows, mNumColumns )
    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );
    SCAI_ASYNCHRONOUS( *syncToken )
    // make all needed data available at loc
    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );
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

    if ( localStorage.getFormat() == Format::ELL )
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
    ContextPtr loc = this->getContextPtr();
    jacobiHalo.getSupportedContext( loc );
    {
        SCAI_CONTEXT_ACCESS( loc )
        WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
        ReadAccess<ValueType> rLocalDiagonal( localDiagonal, loc );
        ReadAccess<IndexType> haloIA( mIA, loc );
        ReadAccess<IndexType> haloJA( mJA, loc );
        ReadAccess<ValueType> haloValues( mValues, loc );
        ReadAccess<ValueType> rOldHaloSolution( haloOldSolution, loc );
        const IndexType numNonEmptyRows = mRowIndexes.size();

        if ( numNonEmptyRows != 0 )
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

    if ( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::asum<ValueType> > asum;
    ContextPtr loc = this->getContextPtr();
    asum.getSupportedContext( loc );
    ReadAccess<ValueType> data( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc );
    return asum[loc]( mValues.size(), data.get(), 1 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    if ( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    static LAMAKernel<blaskernel::BLASKernelTrait::dot<ValueType> > dot;
    ContextPtr loc = this->getContextPtr();
    dot.getSupportedContext( loc );
    ReadAccess<ValueType> data( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc );
    return common::Math::sqrt( dot[loc]( mValues.size(), data.get(), 1, data.get(), 1 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType ELLStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    if ( mNumRows == 0 || mNumValuesPerRow == 0 )
    {
        return static_cast<ValueType>( 0.0 );
    }

    static LAMAKernel<ELLKernelTrait::absMaxVal<ValueType> > absMaxVal;
    ContextPtr loc = this->getContextPtr();
    absMaxVal.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ValueType maxval = absMaxVal[loc]( mNumRows, mNumValuesPerRow, ellIA.get(), ellValues.get() );
    return maxval;
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
    // a and b have to be ELL storages, otherwise create temporaries.
    const ELLStorage<ValueType>* ellA = NULL;
    const ELLStorage<ValueType>* ellB = NULL;
    // Define shared pointers in case we need temporaries
    common::shared_ptr<ELLStorage<ValueType> > tmpA;
    common::shared_ptr<ELLStorage<ValueType> > tmpB;

    if ( a.getFormat() == Format::ELL )
    {
        ellA = dynamic_cast<const ELLStorage<ValueType>*>( &a );
        SCAI_ASSERT_DEBUG( ellA, "could not cast to ELLStorage " << a )
    }
    else
    {
        SCAI_UNSUPPORTED( a << ": will be converted to ELL for matrix add" )
        tmpA = common::shared_ptr<ELLStorage<ValueType> >( new ELLStorage<ValueType>( a ) );
        ellA = tmpA.get();
    }

    if ( b.getFormat() == Format::ELL )
    {
        ellB = dynamic_cast<const ELLStorage<ValueType>*>( &b );
        SCAI_ASSERT_DEBUG( ellB, "could not cast to ELLStorage " << b )
    }
    else
    {
        SCAI_UNSUPPORTED( b << ": will be converted to ELL for matrix add" )
        tmpB = common::shared_ptr<ELLStorage<ValueType> >( new ELLStorage<ValueType>( b ) );
        ellB = tmpB.get();
    }

    ContextPtr loc = this->getContextPtr(); // preferred location for matrix add

    matrixAddMatrixELL( alpha, *ellA, beta, *ellB );
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
    //    common::shared_ptr<ELLStorage<ValueType> > tmpA;
    //    common::shared_ptr<ELLStorage<ValueType> > tmpB;
    common::shared_ptr<ELLStorage<ValueType> > tmpC;

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

    if ( beta != scai::common::constants::ZERO )
    {
        if ( ( c.getFormat() == Format::ELL ) && ( &c != this ) )
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

    if ( beta != scai::common::constants::ZERO )
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
        mNumValuesPerRow = reduce[ loc ]( cIA.get(), mNumRows, 0, utilskernel::binary::MAX );
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
    ContextPtr loc = this->getContextPtr();
    matrixAddSizes.getSupportedContext( loc, reduce, matrixAdd );
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
        mNumValuesPerRow = reduce[loc]( cIA.get(), mNumRows, 0, utilskernel::binary::MAX );
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

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( ELLStorage, SCAI_NUMERIC_TYPES_HOST )

#define ELL_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                 \
    template void ELLStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                              \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                               \
    template void ELLStorage<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;              \
    template void ELLStorage<ValueType>::setRowImpl( const hmemo::HArray<OtherValueType>&, const IndexType,                \
                                                     const utilskernel::binary::BinaryOp );                          \
    template void ELLStorage<ValueType>::getColumnImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;           \
    template void ELLStorage<ValueType>::setColumnImpl( const hmemo::HArray<OtherValueType>&, const IndexType,             \
                                                        const utilskernel::binary::BinaryOp );                       \
    template void ELLStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                          \
    template void ELLStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                          \
    template void ELLStorage<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                                \
    template void ELLStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;                                               \
    template void ELLStorage<ValueType>::setDIADataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );

#define ELL_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, ELL_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

    SCAI_COMMON_LOOP( ELL_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef ELL_STORAGE_INST_LVL2
#undef ELL_STORAGE_INST_LVL1

    } /* end namespace lama */

    } /* end namespace scai */
