/**
 * @file DIAStorage.cpp
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
 * @brief Instantiation for template class DIAStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/DIAStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/DIAKernelTrait.hpp>
#include <scai/sparsekernel/CSRKernelTrait.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/unique_ptr.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

using namespace scai::hmemo;

namespace scai
{

using common::scoped_array;
using common::shared_ptr;

using tasking::SyncToken;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;

using sparsekernel::CSRKernelTrait;
using sparsekernel::DIAKernelTrait;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DIAStorage<ValueType>::logger, "MatrixStorage.DIAStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const IndexType numRows, const IndexType numColumns )

    : CRTPMatrixStorage<DIAStorage<ValueType>, ValueType>( numRows, numColumns ), mNumDiagonals( 0 )
{
    SCAI_LOG_DEBUG( logger, "DIAStorage for matrix " << mNumRows << " x " << mNumColumns << ", no non-zero elements" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage()
    : CRTPMatrixStorage<DIAStorage<ValueType>, ValueType>( 0, 0 ), mNumDiagonals( 0 )
{
    SCAI_LOG_DEBUG( logger, "DIAStorage, matrix is 0 x 0." )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numDiagonals,
    const HArray<IndexType>& offsets,
    const HArray<ValueType>& values )

    : CRTPMatrixStorage<DIAStorage<ValueType>, ValueType>( numRows, numColumns ), mNumDiagonals(
        numDiagonals ), mOffset( offsets ), mValues( values )
{
    // set diagonal property inherited as given
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const DIAStorage<ValueType>& other )

    : CRTPMatrixStorage<DIAStorage<ValueType>, ValueType>( 0, 0 )
{
    // @todo: copy of same storage format should be implemented more efficiently
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "DIAStorage " << mNumRows << " x " << mNumColumns
           << ", #diags = " << mNumDiagonals
           << ", #values = " << mValues.size() << endl;
    ReadAccess<IndexType> offset( mOffset );
    ReadAccess<ValueType> values( mValues );
    stream << "Diagonal offsets:";

    for ( IndexType d = 0; d < mNumDiagonals; d++ )
    {
        stream << " " << offset[d];
    }

    stream << endl;

    for ( IndexType i = 0; i < mNumRows; i++ )
    {
        stream << "Row " << i << " :";

        for ( IndexType ii = 0; ii < mNumDiagonals; ++ii )
        {
            const IndexType j = i + offset[ii];

            if ( !common::Utils::validIndex( j, mNumColumns ) )
            {
                continue;
            }

            stream << " " << j << ":" << values[i + ii * mNumRows];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::clear()
{
    mNumRows = 0;
    mNumColumns = 0;
    mNumDiagonals = 0;
    mOffset.clear();
    mValues.clear();
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format::MatrixStorageFormat DIAStorage<ValueType>::getFormat() const
{
    return Format::DIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setDiagonalImpl( const ValueType value )
{
    SCAI_ASSERT_ERROR( mDiagonalProperty, *this << ": has not diagonal property, cannot set diagonal" )
    IndexType numDiagonalElements = common::Math::min( mNumColumns, mNumRows );
    static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
    // take context of this storage to set
    ContextPtr loc = this->getContextPtr();
    setVal.getSupportedContext( loc );
    {
        // not all values might be changed, so use WriteAccess instead of WriteOnlyAccess
        WriteAccess<ValueType> wValues( mValues, loc );
        ReadAccess<IndexType> rOffset( mOffset, loc );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[loc]( wValues.get(), numDiagonalElements, value, common::binary::COPY );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getRowImpl( HArray<OtherType>& row, const IndexType i ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )

    WriteOnlyAccess<OtherType> wRow( row, mNumColumns );

    const ReadAccess<IndexType> offset( mOffset );
    const ReadAccess<ValueType> values( mValues );

    #pragma omp parallel for

    for ( IndexType j = 0; j < mNumColumns; ++j )
    {
        wRow[j] = static_cast<OtherType>( 0 );
    }

    #pragma omp parallel for

    for ( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( common::Utils::validIndex( j, mNumColumns ) )
        {
            wRow[j] = static_cast<OtherType>( values[diaindex( i, d, mNumRows, mNumDiagonals )] );
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::_HArray& values, const IndexType i ) const
{
    if ( values.getValueType() != this->getValueType() )
    {
        SCAI_LOG_WARN( logger, "use temporary HArray<" << this->getValueType() << ">"
                               << " for getSparseRow with values : HArray<" << values.getValueType() << ">" )
        HArray<ValueType> typedValues;
        getSparseRow( jA, typedValues, i );
        HArrayUtils::assign( values, typedValues );
        return;
    }

    SCAI_REGION( "Storage.DIA.getSparseRow" )

    HArray<ValueType>& typedValues = reinterpret_cast<HArray<ValueType>&>( values );

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )

    WriteOnlyAccess<IndexType> wJA( jA, mNumDiagonals );
    WriteOnlyAccess<ValueType> wValues( typedValues, mNumDiagonals );

    const ReadAccess<IndexType> offset( mOffset );
    const ReadAccess<ValueType> dValues( mValues );

    IndexType n = 0;  // count the real number of entries, might be lt mNumDiagonals 

    for ( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( common::Utils::validIndex( j, mNumColumns ) )
        {
            wJA[n] = j;
            wValues[n] = dValues[diaindex( i, d, mNumRows, mNumDiagonals )];
            ++n;
        }
    }

    wJA.resize( n );
    wValues.resize( n );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getSparseColumn(
    hmemo::HArray<IndexType>& iA,
    hmemo::_HArray& values,
    const IndexType j ) const
{
    SCAI_REGION( "Storage.DIA.getSparseCol" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "col index out of range" )

    ContextPtr loc = Context::getHostPtr();  // only on host here

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        WriteOnlyAccess<IndexType> wRowIndexes( iA, loc, mNumRows );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, mNumRows );

        const ReadAccess<IndexType> offset( mOffset );
        const ReadAccess<ValueType> values( mValues );

        IndexType cnt = 0;

        for ( IndexType d = 0; d < mNumDiagonals; ++d )
        {
            IndexType i = j - offset[d];

            if ( common::Utils::validIndex( i, mNumRows ) )
            {
                wRowIndexes[ cnt ] = i;
                wValuePos  [ cnt ] = diaindex( i, d, mNumRows, mNumDiagonals );
                cnt++;
            }
        }

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    HArrayUtils::gather( values, mValues, valuePos, common::binary::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::setRowImpl( const HArray<OtherType>& row, const IndexType i,
                                        const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), mNumColumns, "row array to small for set" )

    ReadAccess<OtherType> rRow( row );

    ReadAccess<IndexType> offset( mOffset );
    WriteAccess<ValueType> values( mValues );

    #pragma omp parallel for

    for ( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( common::Utils::validIndex( j, mNumColumns ) )
        {
            ValueType val = static_cast<ValueType>( rRow[j] );
            ValueType& loc = values[diaindex( i, d, mNumRows, mNumDiagonals )];

            switch ( op )
            {
                case common::binary::COPY :
                    loc = val;
                    break;
                case common::binary::ADD  :
                    loc += val;
                    break;
                case common::binary::SUB  :
                    loc -= val;
                    break;
                case common::binary::MULT :
                    loc *= val;
                    break;
                case common::binary::DIVIDE :
                    loc /= val;
                    break;
                default:
                    break;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getColumn( _HArray& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.DIA.getDenseCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, mNumRows, colValues, rowIndexes );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::setColumnImpl( const HArray<OtherType>& column, const IndexType j,
        const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), mNumRows, "column array to small for set" )

    ReadAccess<OtherType> rColumn( column );

    ReadAccess<IndexType> offset( mOffset );
    WriteAccess<ValueType> values( mValues );

    #pragma omp parallel for

    for ( IndexType d = 0; d < mNumDiagonals; ++d )
    {
        IndexType i = j - offset[d];

        if ( common::Utils::validIndex( i, mNumRows ) )
        {
            ValueType val = static_cast<ValueType>( rColumn[i] );
            ValueType& loc = values[diaindex( i, d, mNumRows, mNumDiagonals )];

            switch ( op )
            {
                case common::binary::COPY :
                    loc = val;
                    break;
                case common::binary::ADD  :
                    loc += val;
                    break;
                case common::binary::SUB  :
                    loc -= val;
                    break;
                case common::binary::MULT :
                    loc *= val;
                    break;
                case common::binary::DIVIDE :
                    loc /= val;
                    break;
                default:
                    break;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::getDiagonalImpl( HArray<OtherType>& diagonal ) const
{
    static LAMAKernel<UtilKernelTrait::set<OtherType, ValueType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    IndexType numDiagonalElements = common::Math::min( mNumColumns, mNumRows );
    WriteOnlyAccess<OtherType> wDiagonal( diagonal, loc, numDiagonalElements );
    ReadAccess<ValueType> rValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    // Diagonal is first column
    set[ loc ]( wDiagonal.get(), rValues.get(), numDiagonalElements, common::binary::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::setDiagonalImpl( const HArray<OtherType>& diagonal )
{
    IndexType numDiagonalElements = common::Math::min( mNumColumns, mNumRows );
    numDiagonalElements = common::Math::min( numDiagonalElements, diagonal.size() );
    static LAMAKernel<UtilKernelTrait::set<ValueType, OtherType> > set;
    ContextPtr loc = this->getContextPtr();
    set.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<OtherType> rDiagonal( diagonal, loc );
    WriteAccess<ValueType> wValues( mValues, loc );
    set[loc]( wValues.get(), rDiagonal.get(), numDiagonalElements, common::binary::COPY );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scaleImpl( const ValueType value )
{
    HArrayUtils::compute( mValues, mValues, common::binary::MULT, value, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::unary::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherType>
void DIAStorage<ValueType>::scaleImpl( const HArray<OtherType>& diagonal )
{
    {
        ReadAccess<OtherType> rDiagonal( diagonal );
        WriteAccess<ValueType> wValues( mValues );
        ReadAccess<IndexType> rOffset( mOffset );

        for ( IndexType i = 0; i < mNumRows; i++ )
        {
            for ( IndexType ii = 0; ii < mNumDiagonals; ++ii )
            {
                const IndexType j = i + rOffset[ii];

                if ( common::Utils::validIndex( j, mNumColumns ) )
                {
                    wValues[ii * mNumRows + i] *= static_cast<ValueType>( rDiagonal[j] );
                }
            }
        }
    }

    if ( SCAI_LOG_TRACE_ON( logger ) )
    {
        SCAI_LOG_TRACE( logger, "DIA after scale diagonal" )
        print();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
bool DIAStorage<ValueType>::checkDiagonalProperty() const
{
    bool diagonalProperty = true;

    if ( mNumRows != mNumColumns )
    {
        diagonalProperty = false;
    }
    else if ( mNumRows == 0 )
    {
        // zero sized matrix has diagonal property
        diagonalProperty = true;
    }
    else if ( mOffset.size() == 0 )
    {
        // full zero matrix but not zero size -> no diagonal property
        diagonalProperty = false;
    }
    else
    {
        // diagonal property is given if first diagonal is the main one
        ReadAccess<IndexType> offset( mOffset );
        diagonalProperty = offset[0] == 0;
    }

    SCAI_LOG_INFO( logger, *this << ": checkDiagonalProperty -> " << diagonalProperty )
    return diagonalProperty;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::check( const char* /* msg */ ) const
{
    SCAI_ASSERT_EQ_ERROR( mNumDiagonals, mOffset.size(), "size mismatch for array with diagonal offsets" )

    if ( mNumDiagonals == 0 && mNumRows > 0 )
    {
        SCAI_ASSERT_EQUAL_ERROR( false, mDiagonalProperty )
    }

    SCAI_ASSERT_EQUAL_ERROR( mNumDiagonals * mNumRows, mValues.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_DEBUG( logger, "set identity, size = " << size )
    mNumRows = size;
    mNumColumns = size;
    mNumDiagonals = 1; // identity has exactly one diagonal
    {
        static LAMAKernel<UtilKernelTrait::setVal<IndexType> > setVal;
        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );
        WriteOnlyAccess<IndexType> wOffset( mOffset, loc, mNumDiagonals );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[ loc ]( wOffset.get(), 1, 0, common::binary::COPY );
    }
    {
        static LAMAKernel<UtilKernelTrait::setVal<ValueType> > setVal;
        ContextPtr loc = this->getContextPtr();
        setVal.getSupportedContext( loc );
        WriteOnlyAccess<ValueType> values( mValues, loc, mNumRows );
        SCAI_CONTEXT_ACCESS( loc )
        setVal[ loc ]( values.get(), mNumRows, ValueType( 1 ), common::binary::COPY );
    }
    mDiagonalProperty = true;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setUsedDiagonal(
    bool upperDiagonalUsed[],
    bool lowerDiagonalUsed[],
    IndexType i,
    IndexType j )
{
    if ( j >= i )
    {
        bool& flag = upperDiagonalUsed[j - i];

        // set flag only if not already set, improves cache usage

        if ( !flag )
        {
            flag = true; // write only if not true,
        }
    }
    else
    {
        bool& flag = lowerDiagonalUsed[i - j];

        if ( !flag )
        {
            flag = true; // write only if not true,
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename CSRValueType>
void DIAStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<CSRValueType>* values,
    const ContextPtr prefLoc ) const
{
    SCAI_REGION( "Storage.DIA.buildCSR" )
    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;
    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, CSRValueType> > getCSRValues;
    // do it where all routines are avaialble
    ContextPtr loc = prefLoc;
    sizes2offsets.getSupportedContext( loc, getCSRSizes, getCSRValues );
    SCAI_LOG_INFO( logger,
                   "buildTypedCSRData<" << common::getScalarType<CSRValueType>() << ">"
                   << " from DIA<" << common::getScalarType<ValueType>() << "> = " << *this << ", diagonal property = " << mDiagonalProperty )
    ReadAccess<IndexType> diaOffsets( mOffset );
    ReadAccess<ValueType> diaValues( mValues );
    WriteOnlyAccess<IndexType> csrIA( ia, loc, mNumRows + 1 );
    // In contrary to COO and CSR, the DIA format stores also some ZERO values like Dense
    ValueType eps = static_cast<ValueType>( 0.0 );
    getCSRSizes[loc]( csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(),
                      diaValues.get(), eps );

    if ( ja == NULL || values == NULL )
    {
        csrIA.resize( mNumRows );
        return;
    }

    IndexType numValues = sizes2offsets[loc]( csrIA.get(), mNumRows );
    SCAI_LOG_INFO( logger, "CSR: #non-zero values = " << numValues )
    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<CSRValueType> csrValues( *values, loc, numValues );
    getCSRValues[loc]( csrJA.get(), csrValues.get(), csrIA.get(), mDiagonalProperty, mNumRows, mNumColumns,
                       mNumDiagonals, diaOffsets.get(), diaValues.get(), eps );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getFirstColumnIndexes( hmemo::HArray<IndexType>& ) const
{
    COMMON_THROWEXCEPTION( "getFirstColumnIndexes not possible for DENSE format" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValues,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Storage.DIA.setCSR" )

    if ( ia.size() == numRows )
    {
        // offset array required
        HArray<IndexType> offsets;
        IndexType total = _MatrixStorage::sizes2offsets( offsets, ia, prefLoc );
        SCAI_ASSERT_EQUAL( numValues, total, "sizes do not sum to number of values" );
        setCSRDataImpl( numRows, numColumns, numValues, offsets, ja, values, prefLoc );
        return;
    }

    SCAI_ASSERT_EQUAL_DEBUG( numRows + 1, ia.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, ja.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, values.size() )
    static LAMAKernel<CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;
    // prefLoc is ignored, we do it on the Host
    // ToDo: replace Host code with kernels, implement kernels for other devices
    ContextPtr loc = Context::getHostPtr();
    ReadAccess<IndexType> csrIA( ia, loc );
    ReadAccess<IndexType> csrJA( ja, loc );
    ReadAccess<OtherValueType> csrValues( values, loc );
    _MatrixStorage::setDimension( numRows, numColumns );
    SCAI_LOG_DEBUG( logger, "fill DIA sparse matrix " << mNumRows << " x " << mNumColumns << " from csr data" )
    // build a set of all used lower and upper diagonals
    IndexType maxNumDiagonals = common::Math::max( mNumRows, mNumColumns );
    scoped_array<bool> upperDiagonalUsed( new bool[maxNumDiagonals] );
    scoped_array<bool> lowerDiagonalUsed( new bool[maxNumDiagonals] );

    for ( IndexType i = 0; i < maxNumDiagonals; i++ )
    {
        upperDiagonalUsed[i] = false;
        lowerDiagonalUsed[i] = false;
    }

    #pragma omp parallel for 

    for ( IndexType i = 0; i < mNumRows; ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; jj++ )
        {
            IndexType j = csrJA[jj]; // column used
            setUsedDiagonal( upperDiagonalUsed.get(), lowerDiagonalUsed.get(), i, j );
        }
    }

    mDiagonalProperty = hasDiagonalProperty[loc]( numRows, csrIA.get(), csrJA.get() );
    // mDiagonalProperty forces upper diagonal to be the first one
    setOffsets( maxNumDiagonals, upperDiagonalUsed.get(), lowerDiagonalUsed.get() );
    // now we can allocate and set the values
    {
        ReadAccess<IndexType> offset( mOffset );
        WriteOnlyAccess<ValueType> myValues( mValues, mNumDiagonals * mNumRows );
        #pragma omp parallel for 

        for ( IndexType i = 0; i < mNumRows; i++ )
        {
            for ( IndexType d = 0; d < mNumDiagonals; d++ )
            {
                ValueType& addrValue = myValues[diaindex( i, d, mNumRows, mNumDiagonals )];
                // check for j >= 0 and j < mNumColumns not needed here
                addrValue = ValueType( 0 );
                IndexType j = i + offset[d];

                if ( !common::Utils::validIndex( j, mNumColumns ) )
                {
                    continue;
                }

                for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
                {
                    if ( csrJA[jj] == j )
                    {
                        addrValue = static_cast<ValueType>( csrValues[jj] );
                        break;
                    }
                }
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::setDIADataImpl(
    const IndexType /*numRows*/,
    const IndexType /*numColumns*/,
    const IndexType /*numDiagonals*/,
    const HArray<IndexType>& /*offsets*/,
    const HArray<OtherValueType>& /*values*/,
    const ContextPtr /*prefLoc*/ )
{
    COMMON_THROWEXCEPTION( "not yet implemeted" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setOffsets(
    const IndexType maxNumDiagonals,
    const bool upperDiagonalUsed[],
    const bool lowerDiagonalUsed[] )
{
    SCAI_LOG_INFO( logger, "setOffsets, max diagonals = " << maxNumDiagonals << ", " << *this )

    // Help routine to set offsets from used diagonals
    mNumDiagonals = 0;
    IndexType firstIndex = 0;

    if ( mDiagonalProperty )
    {
        firstIndex = 1;
        mNumDiagonals = 1;
    }

    for ( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
    {
        if ( upperDiagonalUsed[i] )
        {
            mNumDiagonals++;
        }

        if ( lowerDiagonalUsed[i] )
        {
            mNumDiagonals++;
        }
    }

    SCAI_LOG_INFO( logger, "storage data requires " << mNumDiagonals << " diagonals a " << mNumRows << " values" )

    WriteOnlyAccess<IndexType> wOffset( mOffset, mNumDiagonals );

    if ( mNumDiagonals > 0 )
    {
        mNumDiagonals = 0;
        firstIndex = 0;

        if ( mDiagonalProperty )
        {
            wOffset[mNumDiagonals++] = 0;
            firstIndex = 1;
        }

        if ( maxNumDiagonals )  // otherwise maxNumDiagonals - 1 is illegal for unsigned
        {
            for ( IndexType i = maxNumDiagonals - 1; i > 0; i-- )
            {
                if ( lowerDiagonalUsed[i] )
                {
                    wOffset[mNumDiagonals++] = -i;
                }
            }
        }

        SCAI_LOG_INFO( logger, "lower diagonals = " << mNumDiagonals )

        for ( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
        {
            if ( upperDiagonalUsed[i] )
            {
                wOffset[mNumDiagonals++] = i;
            }
        }

        SCAI_LOG_INFO( logger, "lower + upper diagonals = " << mNumDiagonals )
    }

    SCAI_ASSERT_EQUAL_DEBUG( mNumDiagonals, wOffset.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::~DIAStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~DIAStorage for matrix " << mNumRows << " x " << mNumColumns << ", # diags = " << mNumDiagonals )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::purge()
{
    mNumColumns = 0;
    mNumRows = 0;
    mNumDiagonals = 0;
    mOffset.purge();
    mValues.purge();
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate DIA sparse matrix of size " << numRows << " x " << numColumns )
    clear();
    mNumRows = numRows;
    mNumColumns = numColumns;
    mDiagonalProperty = checkDiagonalProperty();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DIAStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << mNumRows << " x " << mNumColumns
           << ", nd = " << mNumDiagonals << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )
    ContextPtr prefLoc = this->getContextPtr();
    return utilskernel::HArrayUtils::asum( mValues, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )
    ContextPtr prefLoc = this->getContextPtr();
    ValueType res = HArrayUtils::dotProduct( mValues, mValues, prefLoc );
    return common::Math::sqrt( res );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
typename DIAStorage<ValueType>::StorageAbsType DIAStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )
    static LAMAKernel<DIAKernelTrait::absMaxVal<ValueType> > absMaxVal;
    ContextPtr loc = this->getContextPtr();
    absMaxVal.getSupportedContext( loc );
    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    StorageAbsType maxval = absMaxVal[loc]( mNumRows, mNumColumns, mNumDiagonals, diaOffsets.get(), diaValues.get() );
    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<DIAKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rOffset( mOffset, loc );

    IndexType pos = getValuePos[loc]( i, j, mNumRows, rOffset.get(), mNumDiagonals );

    ValueType val = 0;

    if ( pos != nIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumRows * mNumDiagonals, "illegal value position for ( " << i << ", " << j << " )" );

        val = utilskernel::HArrayUtils::getVal<ValueType>( mValues, pos );
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const common::binary::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, mNumRows, "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, mNumColumns, "column index out of range" )

    SCAI_LOG_TRACE( logger, "set value (" << i << ", " << j << ")" )

    static LAMAKernel<DIAKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rOffset( mOffset, loc );

    IndexType pos = getValuePos[loc]( i, j, mNumRows, rOffset.get(), mNumDiagonals );

    if ( pos == nIndex )
    {
        COMMON_THROWEXCEPTION( "DIA storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mNumRows * mNumDiagonals, "illegal value position for ( " << i << ", " << j << " )" );

    utilskernel::HArrayUtils::setValImpl( mValues, pos, val, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::prefetch( const ContextPtr location ) const
{
    mOffset.prefetch( location );
    mValues.prefetch( location );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DIAStorage<ValueType>::getNumDiagonals() const
{
    return mNumDiagonals;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<IndexType>& DIAStorage<ValueType>::getOffsets() const
{
    return mOffset;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
const HArray<ValueType>& DIAStorage<ValueType>::getValues() const
{
    return mValues;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::wait() const
{
    mOffset.wait();
    mValues.wait();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::swap( _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( getFormat(), other.getFormat(), "swap only for same storage format" )
    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(), "swap only for same value type" )

    // only in debug mode use the more expensive dynamic cast for verification

    SCAI_ASSERT_DEBUG( dynamic_cast<DIAStorage<ValueType>* >( &other ), "illegal storage to swap" )

    swapImpl( reinterpret_cast<DIAStorage<ValueType>& >( other ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::swapImpl( DIAStorage<ValueType>& other )
{
    std::swap( mNumDiagonals, other.mNumDiagonals );
    mOffset.swap( other.mOffset );
    mValues.swap( other.mValues );
    MatrixStorage<ValueType>::swapMS( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t DIAStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( IndexType );
    memoryUsage += sizeof( IndexType ) * mOffset.size();
    memoryUsage += sizeof( ValueType ) * mValues.size();
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.DIA.timesVector" )
    SCAI_LOG_INFO( logger,
                   "Computing z = " << alpha << " * A * x + " << beta << " * y"
                   << ", with A = " << *this << ", x = " << x << ", y = " << y << ", z = " << result )

    if ( alpha == common::constants::ZERO )
    {
        // so we just have result = beta * y, will be done synchronously
        HArrayUtils::compute( result, beta, common::binary::MULT, y, this->getContextPtr() );
        return;
    }

    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )

    if ( beta != common::constants::ZERO )
    {
        SCAI_ASSERT_EQUAL( y.size(), mNumRows, "size mismatch y, beta = " << beta )
    }

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;
    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    // Note: read access to y must appear before write access to result in case of alias
    ReadAccess<ValueType> rX( x, loc );

    if ( beta != common::constants::ZERO )
    {
        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );  // result might be aliased to y
        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                         diaOffsets.get(), diaValues.get() );
    }
    else
    {
        // do not access y at all
        WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );  // result might be aliased to y
        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, NULL, mNumRows, mNumColumns, mNumDiagonals,
                         diaOffsets.get(), diaValues.get() );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::vectorTimesMatrix(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrix, result = " << result << ", alpha = " << alpha << ", x = " << x << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.DIA.VectorTimesMatrix" )
    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumRows )
    ContextPtr loc = this->getContextPtr();

    // Due to DIA format GEVM does not benefit of coupling all in one operation, so split it

    // Step 1: result = beta * y

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumColumns );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        // Note: assignScaled will deal with
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
        HArrayUtils::compute( result, beta, common::binary::MULT, y, loc );
    }

    // Step 2: result = alpha * x * this + 1 * result
    bool async = false;
    SyncToken* token = incGEVM( result, alpha, x, async );
    SCAI_ASSERT( NULL == token, "syncrhonous execution has no token" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_REGION( "Storage.DIA.timesVectorAsync" )
    ContextPtr loc = this->getContextPtr();
    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;
    normalGEMV.getSupportedContext( loc );
    // logging + checks not needed when started as a task
    SCAI_LOG_INFO( logger,
                   "Start z = " << alpha << " * A * x + " << beta << " * y, with A = " << *this
                   << ", x = " << x << ", y = " << y << ", z = " << result << " on " << *loc )
    SCAI_ASSERT_EQUAL_ERROR( x.size(), mNumColumns )
    SCAI_ASSERT_EQUAL_ERROR( y.size(), mNumRows )
    common::unique_ptr<SyncToken> syncToken( loc->getSyncToken() );
    SCAI_ASYNCHRONOUS( *syncToken )
    // all accesses will be pushed to the sync token as LAMA arrays have to be protected up
    // to the end of the computations.
    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, mNumRows );
    SCAI_CONTEXT_ACCESS( loc )
    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), mNumRows, mNumColumns, mNumDiagonals,
                     diaOffsets.get(), diaValues.get() );
    syncToken->pushRoutine( rY.releaseDelayed() );
    syncToken->pushRoutine( wResult.releaseDelayed() );
    syncToken->pushRoutine( rX.releaseDelayed() );
    syncToken->pushRoutine( diaValues.releaseDelayed() );
    syncToken->pushRoutine( diaOffsets.releaseDelayed() );
    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::vectorTimesMatrixAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger,
                   *this << ": vectorTimesMatrixAsync, result = " << result << ", alpha = " << alpha << ", x = " << x
                   << ", beta = " << beta << ", y = " << y )
    SCAI_REGION( "Storage.DIA.vectorTimesMatrixAsync" )
    // Step 1: result = beta * y
    ContextPtr loc = this->getContextPtr();

    if ( beta == common::constants::ZERO )
    {
        result.clear();
        result.resize( mNumColumns );
        HArrayUtils::setScalar( result, ValueType( 0 ), common::binary::COPY, loc );
    }
    else
    {
        // Note: assignScaled will deal with
        SCAI_ASSERT_EQUAL( y.size(), mNumColumns, "size mismatch y, beta = " << beta )
        HArrayUtils::compute( result, beta, common::binary::MULT, y, loc );
    }

    bool async = true;
    // Step 2: result = beta * y
    return incGEVM( result, alpha, x, async );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::incGEVM(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    bool async ) const
{
    SCAI_LOG_INFO( logger, "incGEVM ( async = " << async << " ) , result += " << alpha << " * x * storage" )
    static LAMAKernel<DIAKernelTrait::normalGEVM<ValueType> > normalGEVM;
    ContextPtr loc = this->getContextPtr();
    normalGEVM.getSupportedContext( loc );
    common::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues(  mValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc, mNumColumns );
    // use general kernel, might change
    ValueType beta = 1;
    normalGEVM[loc]( wResult.get(), alpha, rX.get(), beta, wResult.get(),
                     mNumRows, mNumColumns, mNumDiagonals,
                     diaOffsets.get(), diaValues.get() );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( diaOffsets.releaseDelayed() );
        syncToken->pushRoutine( diaValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
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
    static LAMAKernel<DIAKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> diaOffset( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, mNumRows );
    jacobi[loc]( wSolution.get(), mNumColumns, mNumDiagonals, diaOffset.get(), diaValues.get(),
                 rOldSolution.get(), rRhs.get(), omega, mNumRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::copy() const
{
    return new DIAStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::newMatrixStorage() const
{
    common::unique_ptr<DIAStorage<ValueType> > storage( new DIAStorage<ValueType>() );
    storage->setContextPtr( this->getContextPtr() );
    return storage.release();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
_MatrixStorage* DIAStorage<ValueType>::create()
{
    return new DIAStorage<ValueType>();
}

template<typename ValueType>
std::string DIAStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "DIAStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* DIAStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
MatrixStorageCreateKeyType DIAStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DIA, common::getScalarType<ValueType>() );
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DIAStorage, SCAI_NUMERIC_TYPES_HOST )

#define DIA_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                 \
    template void DIAStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                              \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                               \
    template void DIAStorage<ValueType>::getRowImpl( hmemo::HArray<OtherValueType>&, const IndexType ) const;              \
    template void DIAStorage<ValueType>::setRowImpl( const hmemo::HArray<OtherValueType>&, const IndexType,                \
            const common::binary::BinaryOp );                          \
    template void DIAStorage<ValueType>::setColumnImpl( const hmemo::HArray<OtherValueType>&, const IndexType,             \
            const common::binary::BinaryOp );                       \
    template void DIAStorage<ValueType>::getDiagonalImpl( hmemo::HArray<OtherValueType>& ) const;                          \
    template void DIAStorage<ValueType>::setDiagonalImpl( const hmemo::HArray<OtherValueType>& );                          \
    template void DIAStorage<ValueType>::scaleImpl( const hmemo::HArray<OtherValueType>& );                                \
    template void DIAStorage<ValueType>::buildCSR( hmemo::HArray<IndexType>&, hmemo::HArray<IndexType>*,                   \
            hmemo::HArray<OtherValueType>*, const hmemo::ContextPtr ) const;                                               \
    template void DIAStorage<ValueType>::setDIADataImpl( const IndexType, const IndexType, const IndexType,                \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );

#define DIA_STORAGE_INST_LVL1( ValueType )                                                                                  \
    SCAI_COMMON_LOOP_LVL2( ValueType, DIA_STORAGE_INST_LVL2, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( DIA_STORAGE_INST_LVL1, SCAI_NUMERIC_TYPES_HOST )

#undef DIA_STORAGE_INST_LVL2
#undef DIA_STORAGE_INST_LVL1

} /* end namespace lama */

} /* end namespace scai */
