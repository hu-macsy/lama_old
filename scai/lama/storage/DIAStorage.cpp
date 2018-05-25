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
#include <scai/lama/storage/CSRStorage.hpp>

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
#include <scai/common/Constants.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/macros/instantiate.hpp>

#include <memory>

using namespace scai::hmemo;

using std::unique_ptr;

namespace scai
{

using tasking::SyncToken;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::HArrayUtils;

using sparsekernel::CSRKernelTrait;
using sparsekernel::DIAKernelTrait;

using common::BinaryOp;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, DIAStorage<ValueType>::logger, "MatrixStorage.DIAStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( ContextPtr ctx ) :

    MatrixStorage<ValueType>( 0, 0, ctx ),
    mOffset( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "DIAStorage( 0 x 0 ) @ " << *ctx )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( IndexType numRows, IndexType numColumns, ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ),
    mOffset( ctx ),
    mValues( ctx )
{
    SCAI_LOG_DEBUG( logger, "DIAStorage( " << getNumRows() << " x " << getNumColumns() << " @ " << *ctx )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage(
    const IndexType numRows,
    const IndexType numColumns,
    HArray<IndexType> offsets,
    HArray<ValueType> values,
    ContextPtr ctx ) :

    MatrixStorage<ValueType>( numRows, numColumns, ctx ), 
    mOffset( std::move( offsets ) ), 
    mValues( std::move( values  ) )
{
    // set diagonal property inherited as given
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( const DIAStorage<ValueType>& other ) : 

    MatrixStorage<ValueType>( other ),

    mOffset( other.mOffset ),
    mValues( other.mValues )
{
    SCAI_LOG_INFO( logger, "copied DIAStorage other = " << other << ", this = " << *this )
}

/* ------------------------------------------------------------------------------------------------------------------ */

template<typename ValueType>
DIAStorage<ValueType>::DIAStorage( DIAStorage<ValueType>&& other ) noexcept :

    MatrixStorage<ValueType>( std::move( other ) ),

    mOffset( std::move( other.mOffset ) ),
    mValues( std::move( other.mValues ) )
{
    // no further checks as we assume other to be a consistent and valid input DIA storage 
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>& DIAStorage<ValueType>::operator=( const DIAStorage<ValueType>& other )
{
    assignDIA( other );
    return *this;
}

template<typename ValueType>
DIAStorage<ValueType>& DIAStorage<ValueType>::operator=( DIAStorage<ValueType>&& other )
{
    // call of move assignment for base class 

    MatrixStorage<ValueType>::moveImpl( std::move( other ) );

    mOffset = std::move( other.mOffset );
    mValues = std::move( other.mValues );

    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::StorageWrapper<DIAStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    ContextPtr ctx = getContextPtr();   // will force a valid copy in this context

    if ( other.getFormat() == Format::DIA )
    {
        // same format conversion, more efficient solution available

        assignDIA( static_cast<const DIAStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRDataImpl( otherCSR.getNumRows(), otherCSR.getNumColumns(),
                        otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues(), ctx );

        SCAI_LOG_INFO( logger, "assignImpl: other CSR = " << other )
    }
    else
    {
        HArray<IndexType>  csrIA( ctx );
        HArray<IndexType>  csrJA( ctx );
        HArray<ValueType>  csrValues( ctx );     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues, ctx );

        SCAI_LOG_INFO( logger, "assignImpl: other = " << other << " -> tmpCSR: ia = " << csrIA 
                                << ", ja = " << csrJA << ", values = " << csrValues << ", this = " << *this )
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::assignDIA( const DIAStorage<OtherValueType>& other )
{   
    if ( static_cast<const _MatrixStorage*>( &other ) == this )
    {   
        SCAI_LOG_DEBUG( logger, typeName() << ": self assign, skipped, storage = " << other )
        return;
    }
    
    auto ctx = getContextPtr();
    
    // both storage have DIA format, we can just copy the corresponding arrays to the right context
    
    _MatrixStorage::_assign( other );     // assign member variables of base class
    
    HArrayUtils::assign( mOffset, other.getOffsets(), ctx );
    HArrayUtils::assign( mValues, other.getValues(), ctx );
    
    SCAI_LOG_DEBUG( logger, "assignDIA: other = " << other << ", this = " << *this )
}


/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;

    IndexType numDiagonals = mOffset.size();

    stream << "DIAStorage " << getNumRows() << " x " << getNumColumns()
           << ", #diags = " << numDiagonals
           << ", #values = " << mValues.size() << endl;
    ReadAccess<IndexType> offset( mOffset );
    ReadAccess<ValueType> values( mValues );
    stream << "Diagonal offsets:";

    for ( IndexType d = 0; d < numDiagonals; d++ )
    {
        stream << " " << offset[d];
    }

    stream << endl;

    for ( IndexType i = 0; i < getNumRows(); i++ )
    {
        stream << "Row " << i << " :";

        for ( IndexType ii = 0; ii < numDiagonals; ++ii )
        {
            const IndexType j = i + offset[ii];

            if ( !common::Utils::validIndex( j, getNumColumns() ) )
            {
                continue;
            }

            stream << " " << j << ":" << values[i + ii * getNumRows()];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::clear()
{
    _MatrixStorage::setDimension( 0, 0 );

    mOffset.clear();
    mValues.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format DIAStorage<ValueType>::getFormat() const
{
    return Format::DIA;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getRow( HArray<ValueType>& row, const IndexType i ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    WriteOnlyAccess<ValueType> wRow( row, getNumColumns() );

    const ReadAccess<IndexType> offset( mOffset );
    const ReadAccess<ValueType> values( mValues );

    #pragma omp parallel for

    for ( IndexType j = 0; j < getNumColumns(); ++j )
    {
        wRow[j] = ValueType( 0 );
    }

    IndexType numDiagonals = mOffset.size();

    #pragma omp parallel for

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( common::Utils::validIndex( j, getNumColumns() ) )
        {
            wRow[j] = values[diaindex( i, d, getNumRows(), numDiagonals )];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_REGION( "Storage.DIA.getSparseRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    IndexType numDiagonals = mOffset.size();

    // numDiagonals is also maximal number of non-zero entries in a row 

    WriteOnlyAccess<IndexType> wJA( jA, numDiagonals );
    WriteOnlyAccess<ValueType> wValues( values, numDiagonals );

    const ReadAccess<IndexType> offset( mOffset );
    const ReadAccess<ValueType> dValues( mValues );

    IndexType n = 0;  // count the real number of entries, might be lt numDiagonals 

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( common::Utils::validIndex( j, getNumColumns() ) )
        {
            wJA[n] = j;
            wValues[n] = dValues[diaindex( i, d, getNumRows(), numDiagonals )];
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
    hmemo::HArray<ValueType>& values,
    const IndexType j ) const
{
    SCAI_REGION( "Storage.DIA.getSparseCol" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "col index out of range" )

    ContextPtr loc = Context::getHostPtr();  // only on host here

    HArray<IndexType> valuePos;     // positions in the values array

    {
        SCAI_CONTEXT_ACCESS( loc )

        IndexType numDiagonals = mOffset.size();

        WriteOnlyAccess<IndexType> wRowIndexes( iA, loc, getNumRows() );
        WriteOnlyAccess<IndexType> wValuePos( valuePos, loc, getNumRows() );

        const ReadAccess<IndexType> offset( mOffset );
        const ReadAccess<ValueType> values( mValues );

        IndexType cnt = 0;

        for ( IndexType d = 0; d < numDiagonals; ++d )
        {
            IndexType i = j - offset[d];

            if ( common::Utils::validIndex( i, getNumRows() ) )
            {
                wRowIndexes[ cnt ] = i;
                wValuePos  [ cnt ] = diaindex( i, d, getNumRows(), numDiagonals );
                cnt++;
            }
        }

        wRowIndexes.resize( cnt );
        wValuePos.resize( cnt );
    }

    HArrayUtils::gather( values, mValues, valuePos, BinaryOp::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    ReadAccess<ValueType> rRow( row );

    IndexType numDiagonals = mOffset.size();

    ReadAccess<IndexType> offset( mOffset );
    WriteAccess<ValueType> values( mValues );

    #pragma omp parallel for

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        IndexType j = i + offset[d];

        if ( common::Utils::validIndex( j, getNumColumns() ) )
        {
            ValueType val = rRow[j];
            ValueType& loc = values[diaindex( i, d, getNumRows(), numDiagonals )];

            switch ( op )
            {
                case BinaryOp::COPY :
                    loc = val;
                    break;
                case BinaryOp::ADD  :
                    loc += val;
                    break;
                case BinaryOp::SUB  :
                    loc -= val;
                    break;
                case BinaryOp::MULT :
                    loc *= val;
                    break;
                case BinaryOp::DIVIDE :
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
void DIAStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_REGION( "Storage.DIA.getDenseCol" )

    HArray<IndexType> rowIndexes;   // row indexes that have entry for column j
    HArray<ValueType> colValues;    // contains the values of entries belonging to column j

    getSparseColumn( rowIndexes, colValues, j );

    HArrayUtils::buildDenseArray( column, getNumRows(), colValues, rowIndexes, ValueType( 0 ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )
    SCAI_ASSERT_GE_DEBUG( column.size(), getNumRows(), "column array to small for set" )

    ReadAccess<ValueType> rColumn( column );

    ReadAccess<IndexType> offset( mOffset );
    WriteAccess<ValueType> values( mValues );

    IndexType numDiagonals = mOffset.size();

    #pragma omp parallel for

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        IndexType i = j - offset[d];

        if ( common::Utils::validIndex( i, getNumRows() ) )
        {
            ValueType val = rColumn[i];
            ValueType& loc = values[diaindex( i, d, getNumRows(), numDiagonals )];

            switch ( op )
            {
                case BinaryOp::COPY :
                    loc = val;
                    break;
                case BinaryOp::ADD  :
                    loc += val;
                    break;
                case BinaryOp::SUB  :
                    loc -= val;
                    break;
                case BinaryOp::MULT :
                    loc *= val;
                    break;
                case BinaryOp::DIVIDE :
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
void DIAStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    IndexType numDiagonalElements = common::Math::min( getNumColumns(), getNumRows() );

    if ( numDiagonalElements == 0 )
    {
        diagonal.clear();
        return;
    }

    // find the index for the main diagonal in array mOffsets

    IndexType mainIndex = getMainIndex();

    SCAI_LOG_INFO( logger, "getDiagonal, #diagElems = " << numDiagonalElements << ", offset[" << mainIndex << "] == 0" )

    if ( mainIndex == invalidIndex )
    {
        // diagonal not available, so set diagonal 

        HArrayUtils::setSameValue<ValueType>( diagonal, numDiagonalElements, 0 );
    }
    else
    {
        IndexType valuesOffset = mainIndex * getNumRows();
        diagonal.resize( numDiagonalElements );
        HArrayUtils::setArraySection( diagonal, 0, 1, mValues, valuesOffset, 1, numDiagonalElements, BinaryOp::COPY );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    IndexType numDiagonalElements = common::Math::min( getNumColumns(), getNumRows() );

    SCAI_ASSERT_EQ_ERROR( diagonal.size(), numDiagonalElements, "serious mismtach for diagonal" )

    if ( numDiagonalElements == 0 )
    {
        return;
    }

    // find the index for the main diagonal in array mOffsets

    IndexType mainIndex = getMainIndex();

    if ( mainIndex == invalidIndex )
    {
        // diagonal not available, throw an exception

        COMMON_THROWEXCEPTION( "cannot set diagonal in DIAStorage, main diagonal not stored" )
    }
    else
    {
        IndexType valuesOffset = mainIndex * getNumRows();
        HArrayUtils::setArraySection( mValues, valuesOffset, 1, diagonal, 0, 1, numDiagonalElements, BinaryOp::COPY );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
IndexType DIAStorage<ValueType>::getMainIndex() const
{
    auto rOffset = hostReadAccess( mOffset );

    for ( IndexType i = 0; i < mOffset.size(); ++i )
    {
        if ( rOffset[i] == 0 )
        {
            return i;
        }
    }

    return invalidIndex;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setDiagonal( const ValueType value )
{
    IndexType numDiagonalElements = common::Math::min( getNumColumns(), getNumRows() );

    // find the index for the main diagonal in array mOffsets

    IndexType mainIndex = getMainIndex();

    if ( mainIndex == invalidIndex )
    {
        // diagonal not available, throw an exception

        COMMON_THROWEXCEPTION( "cannot set diagonal in DIAStorage, main diagonal not stored" )
    }
    else
    {
        IndexType valuesOffset = mainIndex * getNumRows();
        HArrayUtils::fillArraySection( mValues, valuesOffset, 1, value, numDiagonalElements, BinaryOp::COPY );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scale( const ValueType value )
{
    HArrayUtils::compute( mValues, mValues, BinaryOp::MULT, value, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::conj()
{
    HArrayUtils::unaryOp( mValues, mValues, common::UnaryOp::CONJ, this->getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scaleRows( const HArray<ValueType>& diagonal )
{
    const IndexType numDiagonals = getNumDiagonals();

    {
        ReadAccess<ValueType> rDiagonal( diagonal );
        WriteAccess<ValueType> wValues( mValues );
        ReadAccess<IndexType> rOffset( mOffset );

        for ( IndexType i = 0; i < getNumRows(); i++ )
        {
            for ( IndexType ii = 0; ii < numDiagonals; ++ii )
            {
                const IndexType j = i + rOffset[ii];

                if ( common::Utils::validIndex( j, getNumColumns() ) )
                {
                    wValues[ii * getNumRows() + i] *= rDiagonal[i];
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
void DIAStorage<ValueType>::check( const char* /* msg */ ) const
{
    IndexType numDiagonals = mOffset.size();

    SCAI_ASSERT_EQ_ERROR( numDiagonals * getNumRows(), mValues.size(), 
                          "values array in DIA storage has illegal size" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setIdentity( const IndexType size )
{
    SCAI_LOG_DEBUG( logger, "set identity, size = " << size )

    _MatrixStorage::setDimension( size, size );

    // only main diagonal is available
    mOffset = HArray<IndexType>( { 0 }, getContextPtr() );
  
    HArrayUtils::setSameValue( mValues, size, ValueType( 1 ), getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    _MatrixStorage::setDimension( size, size );

    // only main diagonal is available
    mOffset = HArray<IndexType>( { 0 }, getContextPtr() );
    
    HArrayUtils::setArray( mValues, diagonal, common::BinaryOp::COPY, getContextPtr() );
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

    IndexType numDiagonals = mOffset.size();

    static LAMAKernel<CSRKernelTrait::sizes2offsets> sizes2offsets;
    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;
    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType, CSRValueType> > getCSRValues;
    // do it where all routines are avaialble
    ContextPtr loc = prefLoc;
    sizes2offsets.getSupportedContext( loc, getCSRSizes, getCSRValues );
    SCAI_LOG_INFO( logger,
                   "buildTypedCSRData<" << common::getScalarType<CSRValueType>() << ">"
                   << " from DIA<" << common::getScalarType<ValueType>() << "> = " << *this )
    ReadAccess<IndexType> diaOffsets( mOffset );
    ReadAccess<ValueType> diaValues( mValues );
    WriteOnlyAccess<IndexType> csrIA( ia, loc, getNumRows() + 1 );
    // In contrary to COO and CSR, the DIA format stores also some ZERO values like Dense
    getCSRSizes[loc]( csrIA.get(), getNumRows(), getNumColumns(), numDiagonals, diaOffsets.get(), diaValues.get() );

    if ( ja == NULL || values == NULL )
    {
        csrIA.resize( getNumRows() );
        return;
    }

    IndexType numValues = sizes2offsets[loc]( csrIA.get(), getNumRows() );
    SCAI_LOG_INFO( logger, "CSR: #non-zero values = " << numValues )
    WriteOnlyAccess<IndexType> csrJA( *ja, loc, numValues );
    WriteOnlyAccess<CSRValueType> csrValues( *values, loc, numValues );
    getCSRValues[loc]( csrJA.get(), csrValues.get(), csrIA.get(), getNumRows(), getNumColumns(),
                       numDiagonals, diaOffsets.get(), diaValues.get() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void DIAStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Storage.DIA.setCSR" )

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

    SCAI_ASSERT_EQUAL_DEBUG( numRows + 1, ia.size() )
    SCAI_ASSERT_EQUAL_DEBUG( numValues, values.size() )
    static LAMAKernel<CSRKernelTrait::hasDiagonalProperty> hasDiagonalProperty;
    // prefLoc is ignored, we do it on the Host
    // ToDo: replace Host code with kernels, implement kernels for other devices
    ContextPtr loc = Context::getHostPtr();
    ReadAccess<IndexType> csrIA( ia, loc );
    ReadAccess<IndexType> csrJA( ja, loc );
    ReadAccess<OtherValueType> csrValues( values, loc );
    _MatrixStorage::setDimension( numRows, numColumns );
    SCAI_LOG_DEBUG( logger, "fill DIA sparse matrix " << getNumRows() << " x " << getNumColumns() << " from csr data" )
    // build a set of all used lower and upper diagonals
    IndexType maxNumDiagonals = common::Math::max( getNumRows(), getNumColumns() );
    unique_ptr<bool[]> upperDiagonalUsed( new bool[maxNumDiagonals] );
    unique_ptr<bool[]> lowerDiagonalUsed( new bool[maxNumDiagonals] );

    for ( IndexType i = 0; i < maxNumDiagonals; i++ )
    {
        upperDiagonalUsed[i] = false;
        lowerDiagonalUsed[i] = false;
    }

    #pragma omp parallel for 

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; jj++ )
        {
            IndexType j = csrJA[jj]; // column used
            setUsedDiagonal( upperDiagonalUsed.get(), lowerDiagonalUsed.get(), i, j );
        }
    }

    setOffsets( maxNumDiagonals, upperDiagonalUsed.get(), lowerDiagonalUsed.get() );

    IndexType numDiagonals = mOffset.size();   // available after setOffsets

    // now we can allocate and set the values
    {
        ReadAccess<IndexType> offset( mOffset );
        WriteOnlyAccess<ValueType> myValues( mValues, numDiagonals * getNumRows() );
        #pragma omp parallel for 

        for ( IndexType i = 0; i < getNumRows(); i++ )
        {
            for ( IndexType d = 0; d < numDiagonals; d++ )
            {
                ValueType& addrValue = myValues[diaindex( i, d, getNumRows(), numDiagonals )];
                // check for j >= 0 and j < getNumColumns() not needed here
                addrValue = ValueType( 0 );
                IndexType j = i + offset[d];

                if ( !common::Utils::validIndex( j, getNumColumns() ) )
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

    IndexType numDiagonals = 0;
    IndexType firstIndex = 0;

    for ( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
    {
        if ( upperDiagonalUsed[i] )
        {
            numDiagonals++;
        }

        if ( lowerDiagonalUsed[i] )
        {
            numDiagonals++;
        }
    }

    SCAI_LOG_INFO( logger, "storage data requires " << numDiagonals << " diagonals a " << getNumRows() << " values" )

    WriteOnlyAccess<IndexType> wOffset( mOffset, numDiagonals );

    if ( numDiagonals > 0 )
    {
        numDiagonals = 0;
        firstIndex = 0;

        if ( maxNumDiagonals )  // otherwise maxNumDiagonals - 1 is illegal for unsigned
        {
            for ( IndexType i = maxNumDiagonals - 1; i > 0; i-- )
            {
                if ( lowerDiagonalUsed[i] )
                {
                    wOffset[numDiagonals++] = -i;
                }
            }
        }

        SCAI_LOG_INFO( logger, "lower diagonals = " << numDiagonals )

        for ( IndexType i = firstIndex; i < maxNumDiagonals; i++ )
        {
            if ( upperDiagonalUsed[i] )
            {
                wOffset[numDiagonals++] = i;
            }
        }

        SCAI_LOG_INFO( logger, "lower + upper diagonals = " << numDiagonals )
    }

    SCAI_ASSERT_EQUAL_DEBUG( numDiagonals, wOffset.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>::~DIAStorage()
{
    SCAI_LOG_DEBUG( logger,
                    "~DIAStorage for matrix " << getNumRows() << " x " << getNumColumns() << ", # diags = " << mOffset.size() )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::purge()
{
    _MatrixStorage::setDimension( 0, 0 );

    mOffset.purge();
    mValues.purge();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::allocate( IndexType numRows, IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate DIA sparse matrix of size " << numRows << " x " << numColumns )

    clear();

    _MatrixStorage::setDimension( numRows, numColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "DIAStorage<" << common::getScalarType<ValueType>()
           << ">( size = " << getNumRows() << " x " << getNumColumns()
           << ", nd = " << getNumDiagonals() << " )";
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DIAStorage<ValueType>::l1Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l1Norm()" )
    ContextPtr prefLoc = this->getContextPtr();
    return utilskernel::HArrayUtils::l1Norm( mValues, prefLoc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DIAStorage<ValueType>::l2Norm() const
{
    SCAI_LOG_INFO( logger, *this << ": l2Norm()" )
    ContextPtr prefLoc = this->getContextPtr();
    ValueType res = HArrayUtils::dotProduct( mValues, mValues, prefLoc );
    return common::Math::sqrt( res );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DIAStorage<ValueType>::maxNorm() const
{
    SCAI_LOG_INFO( logger, *this << ": maxNorm()" )

    IndexType numDiagonals = getNumDiagonals();

    static LAMAKernel<DIAKernelTrait::absMaxVal<ValueType> > absMaxVal;
    ContextPtr loc = this->getContextPtr();
    absMaxVal.getSupportedContext( loc );

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );

    SCAI_CONTEXT_ACCESS( loc )
    RealType<ValueType> maxval = absMaxVal[loc]( getNumRows(), getNumColumns(), numDiagonals, diaOffsets.get(), diaValues.get() );
    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "get value (" << i << ", " << j << ")" )

    static LAMAKernel<DIAKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rOffset( mOffset, loc );

    IndexType numDiagonals = mOffset.size();
    IndexType pos = getValuePos[loc]( i, j, getNumRows(), rOffset.get(), numDiagonals );

    ValueType val = 0;

    if ( pos != invalidIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, getNumRows() * numDiagonals, "illegal value position for ( " << i << ", " << j << " )" );

        val = utilskernel::HArrayUtils::getVal<ValueType>( mValues, pos );
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setValue( const IndexType i,
                                      const IndexType j,
                                      const ValueType val,
                                      const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    SCAI_LOG_TRACE( logger, "set value (" << i << ", " << j << ")" )

    static LAMAKernel<DIAKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = this->getContextPtr();
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rOffset( mOffset, loc );

    IndexType numDiagonals = mOffset.size();
    IndexType pos = getValuePos[loc]( i, j, getNumRows(), rOffset.get(), numDiagonals );

    if ( pos == invalidIndex )
    {
        COMMON_THROWEXCEPTION( "DIA storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, getNumRows() * numDiagonals, "illegal value position for ( " << i << ", " << j << " )" );

    utilskernel::HArrayUtils::setVal( mValues, pos, val, op );
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
    return mOffset.size();
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
void DIAStorage<ValueType>::swap( DIAStorage<ValueType>& other )
{
    // swap base class

    MatrixStorage<ValueType>::swap( other );

    // swap my member variables

    mOffset.swap( other.mOffset );
    mValues.swap( other.mValues );
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
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    bool async = false; // synchronously execution, no SyncToken required
    SyncToken* token = gemv( result, alpha, x, beta, y, op, async );
    SCAI_ASSERT( token == NULL, "There should be no sync token for synchronous execution" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op,
    bool async ) const
{
    SCAI_REGION( "Storage.DIA.gemv" )

    const IndexType nSource = common::isTranspose( op ) ? getNumRows() : getNumColumns();
    const IndexType nTarget = common::isTranspose( op ) ? getNumColumns() : getNumRows();

    IndexType numDiagonals = mOffset.size();

    SCAI_LOG_INFO( logger,
                   "gemv<" << getValueType() << "> ( op = " << op << ", async = " << async
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

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = this->getContextPtr();
    normalGEMV.getSupportedContext( loc );

    unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> diaOffsets( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    if ( beta != common::Constants::ZERO )
    {
        SCAI_ASSERT_EQ_ERROR( y.size(), nTarget, "y has illegal size" )

        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );  // result might be aliased to y

        SCAI_LOG_INFO( logger, "call kernel normalGEMV( beta = " << beta << ", y[" << nTarget << "] = " << rY.get() << " ) on " << *loc )

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), 
                         getNumRows(), getNumColumns(), numDiagonals,
                         diaOffsets.get(), diaValues.get(), op );

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

        SCAI_LOG_INFO( logger, "call kernel normalGEMV( beta is 0 ) on " << *loc )

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL, 
                         getNumRows(), getNumColumns(), numDiagonals,
                         diaOffsets.get(), diaValues.get(), op );

        if ( async )
        {
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }

    if ( async )
    {
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( diaValues.releaseDelayed() );
        syncToken->pushRoutine( diaOffsets.releaseDelayed() );
    }

    return syncToken.release();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::matrixTimesVectorAsync(
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
void DIAStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), oldSolution.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), rhs.size() )
    SCAI_ASSERT_EQUAL_DEBUG( getNumRows(), getNumColumns() )
    // matrix must be square
    static LAMAKernel<DIAKernelTrait::jacobi<ValueType> > jacobi;
    ContextPtr loc = this->getContextPtr();
    jacobi.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> diaOffset( mOffset, loc );
    ReadAccess<ValueType> diaValues( mValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, getNumRows() );
    jacobi[loc]( wSolution.get(), getNumColumns(), getNumDiagonals(), diaOffset.get(), diaValues.get(),
                 rOldSolution.get(), rRhs.get(), omega, getNumRows() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::copy() const
{
    return new DIAStorage<ValueType>( *this );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
DIAStorage<ValueType>* DIAStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    std::unique_ptr<DIAStorage<ValueType> > storage( new DIAStorage<ValueType>( getContextPtr() ) );
    storage->allocate( numRows, numColumns );
    return storage.release();
}

/* ========================================================================= */
/*  Static fatory methods and related virtual methods                        */
/* ========================================================================= */

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

template<typename ValueType>
const char* DIAStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType DIAStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::DIA, common::getScalarType<ValueType>() );
}

template<typename ValueType>
MatrixStorageCreateKeyType DIAStorage<ValueType>::getCreateValue() const
{
    return createValue();
}

template<typename ValueType>
_MatrixStorage* DIAStorage<ValueType>::create()
{
    return new DIAStorage<ValueType>();
}

/* ========================================================================= */
/*       Template specializations and instantiations                         */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( DIAStorage, SCAI_NUMERIC_TYPES_HOST )

#define DIA_STORAGE_INST_LVL2( ValueType, OtherValueType )                                                                 \
    template void DIAStorage<ValueType>::setCSRDataImpl( const IndexType, const IndexType,                                 \
            const hmemo::HArray<IndexType>&, const hmemo::HArray<IndexType>&,                                              \
            const hmemo::HArray<OtherValueType>&, const hmemo::ContextPtr );                                               \
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
