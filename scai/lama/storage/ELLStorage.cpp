/**
 * @file ELLStorage.cpp
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
 * @brief Instantitions for template class ELLStorage.
 * @author Lauretta Schubert
 * @date 25.05.2011
 */

// hpp
#include <scai/lama/storage/ELLStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/ELLKernelTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/SparseKernelTrait.hpp>
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
using utilskernel::SparseKernelTrait;
using utilskernel::HArrayUtils;

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
ELLStorage<ValueType>::ELLStorage( ContextPtr ctx ) : 

    MatrixStorage<ValueType>( 0, 0, ctx ), 
    mNumValuesPerRow( 0 ),
    mIA( ctx ),
    mJA( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "ELLStorage, default constructor for zero matrix." )

    _MatrixStorage::resetDiagonalProperty();
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

    _MatrixStorage::resetDiagonalProperty();
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

    _MatrixStorage::resetDiagonalProperty();

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

    _MatrixStorage::resetDiagonalProperty();

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
    mDiagonalProperty = checkDiagonalProperty();
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

    mDiagonalProperty = true;

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

    mDiagonalProperty = true; // obviously given for identity matrix

    // Note: we do not build row indexes, no row is empty

    SCAI_LOG_INFO( logger, *this << ": diagonal matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool ELLStorage<ValueType>::checkDiagonalProperty() const
{
    SCAI_LOG_INFO( logger, "checkDiagonalProperty" )

    IndexType numDiagonals = common::Math::min( getNumRows(), getNumColumns() );

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

    _MatrixStorage::setDimension( 0, 0 );

    mNumValuesPerRow = 0;

    mIA.clear();
    mJA.clear();
    mValues.clear();

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

    check( "ELLStorage::buildCSR" );

    // step 1 : compute IA offsets
    IndexType numValues = 0;
    ia.clear();
    ia.reserve( context, getNumRows() + 1 );  // reserve one more entry
    HArrayUtils::assign( ia, mIA, context );

    if ( ja == NULL || values == NULL )
    {
        return;
    }

    numValues = HArrayUtils::scan1( ia, context );
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
    getCSRValues[loc]( csrJA.get(), csrValues.get(), csrIA.get(), getNumRows(), mNumValuesPerRow,
                       ellSizes.get(), ellJA.get(), ellValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void ELLStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    const ContextPtr context )
{
    SCAI_REGION( "Storage.ELL.setCSR" )

    IndexType numValues = ja.size();

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
    std::unique_ptr<HArray<IndexType> > tmpOffsets;
    const HArray<IndexType>* offsets = &ia;

    if ( ia.size() == numRows + 1 )
    {
        ContextPtr loc = context;
        static LAMAKernel<CSRKernelTrait::offsets2sizes > offsets2sizes;
        offsets2sizes.getSupportedContext( loc );
        ReadAccess<IndexType> csrIA( ia, loc );
        WriteOnlyAccess<IndexType> ellSizes( mIA, loc, getNumRows() );
        SCAI_CONTEXT_ACCESS( loc )
        offsets2sizes[ loc ]( ellSizes.get(), csrIA.get(), getNumRows() );
    }
    else if ( ia.size() == numRows )
    {
        HArrayUtils::assign( mIA, ia, context );
        // as the offset array is also needed
        tmpOffsets.reset( ia.copy() );
        IndexType total = HArrayUtils::scan1( *tmpOffsets, context );
        SCAI_ASSERT_EQUAL( total, numValues, "sizes do not sum up correctly" )
        offsets = tmpOffsets.get();
    }
    else
    {
        COMMON_THROWEXCEPTION( "ia array has illegal size " << ia.size() << " for #rows = " << numRows )
    }

    // determine the maximal number of non-zero in one row
    mNumValuesPerRow = HArrayUtils::reduce( mIA, common::BinaryOp::MAX, context );
    SCAI_LOG_DEBUG( logger, "setCSRData, #values/row = " << mNumValuesPerRow )
    //  Now we know the size of the ja and values arrays for the ELL format
    const IndexType dataSize = mNumValuesPerRow * getNumRows();

    if ( getNumRows() > 200 && mNumValuesPerRow > 0 )
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
                           getNumRows(), mNumValuesPerRow,
                           csrIA.get(), csrJA.get(), csrValues.get() );
        SCAI_LOG_DEBUG( logger, " size = " << ellJA.size() )
        IndexType numDiagonals = std::min( getNumRows(), getNumColumns() );

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
            mDiagonalProperty = hasDiagonalProperty[loc]( numDiagonals, ellJA.get() );
        }
    }

    if ( numRows == numColumns && !mDiagonalProperty )
    {
        SCAI_LOG_INFO( logger, *this << ": square matrix has not diagonal property" )
    }

    buildRowIndexes( loc );

    SCAI_LOG_INFO( logger, "ELL: set CSR data done, this = " << *this )
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
    this->resetDiagonalProperty();
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
    SCAI_LOG_INFO( logger, "setDiagonalImpl # value = " << value )
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    ContextPtr loc = this->getContextPtr();
    setVal.getSupportedContext( loc );
    IndexType numDiagonalElements = std::min( getNumColumns(), getNumRows() );
    SCAI_CONTEXT_ACCESS( loc )
    WriteAccess<ValueType> wValues( mValues, loc );
    setVal[ loc ]( wValues.get(), numDiagonalElements, value, common::BinaryOp::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    SCAI_ASSERT_ERROR( hasDiagonalProperty(), "cannot set diagonal for CSR, no diagonal property" )

    const IndexType numDiagonalElements = std::min( getNumColumns(), getNumRows() );

    SCAI_LOG_INFO( logger, "setDiagonalV # diagonal = " << diagonal )
    static LAMAKernel<UtilKernelTrait::set<ValueType, ValueType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    // ELL format with diagonal property: diagonal is just the first column in mValues
    set[ loc ]( wValues.get(), rDiagonal.get(), numDiagonalElements, common::BinaryOp::COPY );
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
void ELLStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const
{   
    SCAI_REGION( "Storage.ELL.getSparseCol" )
    
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )
    
    static LAMAKernel<ELLKernelTrait::getValuePosCol> getValuePosCol;

    ContextPtr loc = this->getContextPtr();

    getValuePosCol.getSupportedContext( loc );

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( iA, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, getNumRows() );

        ReadAccess<IndexType> rIA( mIA, loc );
        ReadAccess<IndexType> rJA( mJA, loc );

        IndexType cnt = getValuePosCol[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), getNumRows(), rJA.get(), mNumValuesPerRow );

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    // column_values = mValues[ pos ];

    HArrayUtils::gather( values, mValues, valuePos, common::BinaryOp::COPY, loc );
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
    SCAI_LOG_INFO( logger, "getDiagonal # diagonal = " << diagonal )
    IndexType numDiagonalElements = common::Math::min( getNumColumns(), getNumRows() );
    static LAMAKernel<UtilKernelTrait::set<ValueType, ValueType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    WriteOnlyAccess<ValueType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );
    // ELL format with diagonal property: diagonal is just the first column in mValues
    SCAI_CONTEXT_ACCESS( loc )
    set[loc]( wDiagonal.get(), rValues.get(), numDiagonalElements, common::BinaryOp::COPY );
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
    SCAI_LOG_INFO( logger, "scaleRows # values = " << values )
    static LAMAKernel<ELLKernelTrait::scaleRows<ValueType> > ellScaleRows;
    ContextPtr loc = this->getContextPtr();
    ellScaleRows.getSupportedContext( loc );
    ReadAccess<ValueType> rValues( values, loc );
    ReadAccess<IndexType> rIa( mIA, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    ellScaleRows[loc]( wValues.get(), getNumRows(), mNumValuesPerRow, rIa.get(), rValues.get() );
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

    _MatrixStorage::resetDiagonalProperty();
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

    static LAMAKernel<ELLKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( mIA, loc );
    ReadAccess<IndexType> rJa( mJA, loc );

    IndexType pos = getValuePos[loc]( i, j, getNumRows(), mNumValuesPerRow, rIa.get(), rJa.get() );

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

    static LAMAKernel<ELLKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );

    IndexType pos = invalidIndex;

    {
        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> rIa( mIA, loc );
        ReadAccess<IndexType> rJa( mJA, loc );

        pos = getValuePos[loc]( i, j, getNumRows(), mNumValuesPerRow, rIa.get(), rJa.get() );

    }

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
void ELLStorage<ValueType>::buildRowIndexes( const ContextPtr context )
{
    SCAI_LOG_INFO( logger, "buildRowIndexes # loc = " << context )
    mRowIndexes.clear();

    if ( getNumRows() == 0 )
    {
        return;
    }

    // Note: compress functionality in HArrayUtils available but we
    // reimplement it here in the same way as compress is optionally done
    // depending on the threshold value

    // Get function pointers for needed kernel routines

    static LAMAKernel<SparseKernelTrait::countNonZeros<IndexType> > countNonZeros;
    static LAMAKernel<SparseKernelTrait::compress<IndexType, IndexType> > compress;

    // choose location where both routines are available

    ContextPtr loc = context;
    countNonZeros.getSupportedContext( loc, compress );

    ReadAccess<IndexType> ellIA( mIA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // count the number of non-zero rows to have a good value for allocation of rowIndexes

    IndexType zero = 0;   // sparse storage uses always the real 0
    IndexType eps  = 0;   // no tolerances used here

    IndexType nonZeroRows = countNonZeros[loc]( ellIA.get(), getNumRows(), zero, eps );

    float usage = float( nonZeroRows ) / float( getNumRows() );

    if ( usage >= mCompressThreshold )
    {
        SCAI_LOG_INFO( logger,
                       "ELLStorage: do not build row indexes, usage = " << usage << " >= " << mCompressThreshold << " ( threshold )" )
        return;
    }

    WriteOnlyAccess<IndexType> rowIndexes( mRowIndexes, loc, nonZeroRows );

    IndexType cnt = compress[loc]( NULL, rowIndexes.get(), ellIA.get(), getNumRows(), zero, eps );

    SCAI_ASSERT_EQ_ERROR( cnt, nonZeroRows, "serious mismatch" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::compress( const RealType<ValueType> eps, const bool keepDiagonal )
{
    SCAI_LOG_INFO( logger, "compress: eps = " << eps )

    ContextPtr loc = this->getContextPtr();
    static LAMAKernel<ELLKernelTrait::compressIA<ValueType> > compressIA;
    static LAMAKernel<UtilKernelTrait::reduce<IndexType> > reduce;
    compressIA.getSupportedContext( loc, reduce );

    IndexType newNumValuesPerRow = invalidIndex;

    HArray<IndexType> newIAArray;
    {
        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> IA( mIA, loc );
        ReadAccess<IndexType> JA( mJA, loc );
        ReadAccess<ValueType> values( mValues, loc );
        // 1. Step: Check for 0 elements and write new IA array
        WriteOnlyAccess<IndexType> newIA( newIAArray, loc, getNumRows() );
        compressIA[loc]( newIA.get(), IA.get(), JA.get(), values.get(), getNumRows(), mNumValuesPerRow, eps, keepDiagonal );
        // 2. Step: compute length of longest row
        newNumValuesPerRow = reduce[ loc ]( newIA.get(), getNumRows(), 0, common::BinaryOp::MAX );
    }

    // Do further steps, if new array could be smaller
    if ( newNumValuesPerRow < mNumValuesPerRow )
    {
        static LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;
        compressValues.getSupportedContext( loc );

        SCAI_CONTEXT_ACCESS( loc )

        // 3. Step: Allocate new JA and Values array
        HArray<ValueType> newValuesArray;
        HArray<IndexType> newJAArray;

        {
            ReadAccess<IndexType> IA( mIA, loc );
            ReadAccess<IndexType> JA( mJA, loc );
            ReadAccess<ValueType> values( mValues, loc );
            WriteOnlyAccess<ValueType> newValues( newValuesArray, loc, getNumRows() * newNumValuesPerRow );
            WriteOnlyAccess<IndexType> newJA( newJAArray, loc, getNumRows() * newNumValuesPerRow );
            // 4. Step: Compute new JA and Values array
            compressValues[loc]( newJA.get(), newValues.get(), newNumValuesPerRow,
                                 IA.get(), JA.get(), values.get(), getNumRows(), mNumValuesPerRow, eps, keepDiagonal );
        }

        mIA = std::move( newIAArray );
        mJA = std::move( newJAArray );
        mValues = std::move( newValuesArray );
        mNumValuesPerRow = newNumValuesPerRow;
    }
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

    if ( alpha == common::Constants::ZERO || ( mNumValuesPerRow == 0 ) )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::compute( result, beta, common::BinaryOp::MULT, y, this->getContextPtr() );

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

#if !defined( SCAI_ASSERT_LEVEL_OFF )
    IndexType nSource = common::isTranspose( op ) ? getNumRows() : getNumColumns();
    IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();
#endif

    SCAI_ASSERT_EQ_ERROR( x.size(), nSource, "gemv: A * x, x has illegal size" )

    if ( beta == common::Constants::ZERO )
    {
        // take version that does not access y at all (can be undefined or aliased to result)
        return normalGEMV( result, alpha, x, op, async );
    }

    // y is relevant, so it must have correct size

    SCAI_ASSERT_EQ_ERROR( y.size(), nTarget, "gemv: A * x + y, y has illegal size" )

    if ( &result == &y && ( beta == common::Constants::ONE ) && ( mRowIndexes.size() > 0 ) )
    {
        // y += A * x,  where only some rows in A are filled, uses more efficient routine
        return sparseGEMV( result, alpha, x, op, async );
    }
    else
    {
        return normalGEMV( result, alpha, x, beta, y, op, async );
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
    const common::MatrixOp op,
    bool async ) const
{
    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    if ( loc != this->getContextPtr() )
    {
        SCAI_LOG_INFO( logger, "normalGEMV not on " << *this->getContextPtr() << " but on " << *loc )
    }
    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    const IndexType nResult = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    SCAI_CONTEXT_ACCESS( loc )
    SCAI_ASYNCHRONOUS( syncToken.get() )
    // Note: alias &result == &y possible
    //       ReadAccess on y before WriteOnlyAccess on result guarantees valid data
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, nResult );
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), 
                     getNumRows(), getNumColumns(), mNumValuesPerRow,
                     ellIA.get(), ellJA.get(), ellValues.get(), op );

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
    const common::MatrixOp op,
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

    const IndexType nResult = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    SCAI_ASYNCHRONOUS( syncToken.get() )
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, nResult );
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), 0, NULL, 
                     getNumRows(), getNumColumns(), mNumValuesPerRow,
                     ellIA.get(), ellJA.get(), ellValues.get(), op );

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
    const common::MatrixOp op,
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
    sparseGEMV[loc]( wResult.get(), alpha, rX.get(), getNumRows(), mNumValuesPerRow, numNonZeroRows,
                     rRowIndexes.get(), ellIA.get(), ellJA.get(), ellValues.get(), op );

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

    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), getNumColumns() )
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
    WriteOnlyAccess<ValueType> wSolution( solution, loc, getNumRows() );
    jacobi[loc] ( wSolution.get(), getNumRows(), mNumValuesPerRow, ellSizes.get(), ellJA.get(), ellValues.get(),
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

    if ( loc->getType() == common::ContextType::Host )
    {
        // used later in OpenMP to generate a TaskSyncToken
        void ( ELLStorage::*jb )(
            HArray<ValueType>&,
            const HArray<ValueType>&,
            const HArray<ValueType>&,
            const ValueType omega ) const
        = &ELLStorage<ValueType>::jacobiIterate;
        using std::bind;
        using std::cref;
        using std::ref;
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
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), getNumColumns() )
    std::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );
    SCAI_ASYNCHRONOUS( *syncToken )
    // make all needed data available at loc
    ReadAccess<IndexType> ellSizes( mIA, loc );
    ReadAccess<IndexType> ellJA( mJA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, getNumRows() );
    SCAI_CONTEXT_ACCESS( loc )
    jacobi[loc]( wSolution.get(), getNumRows(), mNumValuesPerRow, ellSizes.get(), ellJA.get(), ellValues.get(),
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
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& haloOldSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.ELL.jacobiIterateHalo" )
    SCAI_LOG_INFO( logger, "HOST: Jacobi iteration on halo matrix data." )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), localSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumColumns(), haloOldSolution.size() )
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
            jacobiHalo[loc]( wSolution.get(), getNumRows(), rLocalDiagonal.get(), mNumValuesPerRow, haloIA.get(), haloJA.get(),
                             haloValues.get(), haloRowIndexes.get(), numNonEmptyRows, rOldHaloSolution.get(), omega );
        }
        else
        {
            // no row indexes available, computation is done over all rows
            const IndexType numNonEmptyRows = getNumRows();
            jacobiHalo[loc]( wSolution.get(), getNumRows(), rLocalDiagonal.get(), mNumValuesPerRow, haloIA.get(), haloJA.get(),
                             haloValues.get(), NULL, numNonEmptyRows, rOldHaloSolution.get(), omega );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void ELLStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{   
    halo.halo2Global( mJA );
    _MatrixStorage::setDimension( getNumRows(), globalNumColumns );
    _MatrixStorage::resetDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> ELLStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )

    if ( getNumRows() == 0 || mNumValuesPerRow == 0 )
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
RealType<ValueType> ELLStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )

    if ( getNumRows() == 0 || mNumValuesPerRow == 0 )
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
RealType<ValueType> ELLStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    if ( getNumRows() == 0 || mNumValuesPerRow == 0 )
    {
        return RealType<ValueType>( 0 );
    }

    static LAMAKernel<ELLKernelTrait::absMaxVal<ValueType> > absMaxVal;
    ContextPtr loc = this->getContextPtr();
    absMaxVal.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ellIA( mIA, loc );
    ReadAccess<ValueType> ellValues( mValues, loc );
    ValueType maxval = absMaxVal[loc]( getNumRows(), mNumValuesPerRow, ellIA.get(), ellValues.get() );
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
    mDiagonalProperty = ( getNumRows() == getNumColumns() );
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
    //mDiagonalProperty = ( getNumRows() == getNumColumns() );
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

#define ELL_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                \
    template void ELLStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                             \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                              \
    template void ELLStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,  \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;              

#define ELL_STORAGE_INST_LVL1( ValueType )                                                                \
    SCAI_COMMON_LOOP_LVL2( ValueType, ELL_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( ELL_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef ELL_STORAGE_INST_LVL2
#undef ELL_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
