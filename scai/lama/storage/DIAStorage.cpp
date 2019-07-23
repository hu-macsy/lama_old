/**
 * @file DIAStorage.cpp
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
 * @brief Instantiation for template class DIAStorage.
 * @author Thomas Brandes
 * @date 04.06.2011
 */

// hpp
#include <scai/lama/storage/DIAStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/sparsekernel/DIAUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/TransferUtils.hpp>
#include <scai/utilskernel/freeFunction.hpp>

#include <scai/hmemo/ContextAccess.hpp>

#include <scai/tasking/NoSyncToken.hpp>

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

using utilskernel::HArrayUtils;
using utilskernel::TransferUtils;

using sparsekernel::CSRUtils;
using sparsekernel::DIAUtils;

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
    // ToDo: there might be some checks, mainly for offsets: legal values, no doubles
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
    if ( other.getFormat() == Format::DIA )
    {
        // same format conversion, more efficient solution available

        assignDIA( static_cast<const DIAStorage<OtherValueType> & >( other ) );
    }
    else if ( other.getFormat() == Format::CSR )
    {
        const auto otherCSR = static_cast<const CSRStorage<OtherValueType> & >( other );

        setCSRData( otherCSR.getNumRows(), otherCSR.getNumColumns(),
                    otherCSR.getIA(), otherCSR.getJA(), otherCSR.getValues() );

        SCAI_LOG_INFO( logger, "assignImpl: other CSR = " << other )
    }
    else
    {
        HArray<IndexType>  csrIA;
        HArray<IndexType>  csrJA;
        HArray<ValueType>  csrValues;     // might also be OtherValueType, depending on size

        other.buildCSRData( csrIA, csrJA, csrValues );

        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues );

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
void DIAStorage<ValueType>::getSparseRow( 
    hmemo::HArray<IndexType>& jA, 
    hmemo::HArray<ValueType>& values, 
    const IndexType i ) const
{
    SCAI_REGION( "Storage.DIA.getSparseRow" )

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    ContextPtr loc = Context::getHostPtr();  // only on host here

    HArray<IndexType> positions;     // positions in the values array

    DIAUtils::getRowPositions( jA, positions, i, 
                               getNumRows(), getNumColumns(), mOffset, getContextPtr() );

    HArrayUtils::gather( values, mValues, positions, BinaryOp::COPY, loc );
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

    HArray<IndexType> positions;     // positions in the values array

    DIAUtils::getColPositions( iA, positions, j, 
                               getNumRows(), getNumColumns(), mOffset, getContextPtr() );

    HArrayUtils::gather( values, mValues, positions, BinaryOp::COPY, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    HArray<IndexType> ia;          // available indexes of the row
    HArray<IndexType> positions;   // positions in the values array

    DIAUtils::getRowPositions( ia, positions, i,
                               getNumRows(), getNumColumns(), mOffset, getContextPtr() );

    SCAI_ASSERT_EQ_DEBUG( ia.size(), positions.size(), "illegal result arrays for row positions." )

    HArray<ValueType> sparseRow;  // available entries of the row
    
    // for each k :  values[positions[k]] op= row[ia[k]] 

    if ( op == BinaryOp::COPY )
    {
        TransferUtils::copy( mValues, positions, row, ia );
    }
    else
    {
        HArray<ValueType> sparseRow;  // gather entries of row that can be set

        HArrayUtils::gather( sparseRow, row, ia, BinaryOp::COPY, getContextPtr() );
        HArrayUtils::scatter( mValues, positions, true, sparseRow, op, getContextPtr() );
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

    HArray<IndexType> ia;          // available indexes of the columns
    HArray<IndexType> positions;   // positions in the values array

    DIAUtils::getColPositions( ia, positions, j,
                               getNumRows(), getNumColumns(), mOffset, getContextPtr() );

    // for each k :  values[positions[k]] op= column[ia[k]] 

    if ( op == BinaryOp::COPY )
    {
        TransferUtils::copy( mValues, positions, column, ia );
    }
    else
    {
        HArray<ValueType> sparseCol;  // gather entries of column that can be set

        HArrayUtils::gather( sparseCol, column, ia, BinaryOp::COPY, getContextPtr() );
        HArrayUtils::scatter( mValues, positions, true, sparseCol, op, getContextPtr() );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    const IndexType numDiagonalElements = this->getDiagonalSize();

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
    const IndexType numDiagonalElements = this->getDiagonalSize();

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
    const IndexType numDiagonalElements = this->getDiagonalSize();

    if ( numDiagonalElements == 0 )
    {
        return;     // there will be no diagonals available at all
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
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::scaleColumns( const HArray<ValueType>& diagonal )
{
    DIAUtils::setColumns( mValues, getNumRows(), getNumColumns(), mOffset,
                          diagonal, common::BinaryOp::MULT, getContextPtr() );
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
void DIAStorage<ValueType>::buildCSRSizes( hmemo::HArray<IndexType>& ia ) const
{
    DIAUtils::getCSRSizes( ia, getNumRows(), getNumColumns(), mOffset, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::buildCSRData(
    hmemo::HArray<IndexType>& csrIA,
    hmemo::HArray<IndexType>& csrJA,
    hmemo::_HArray& csrValues ) const
{
    if ( csrValues.getValueType() == getValueType() )
    {
        HArray<ValueType>& typedCSRValues = static_cast<HArray<ValueType>&>( csrValues );

        DIAUtils::convertDIA2CSR( csrIA, csrJA, typedCSRValues,
                                  getNumRows(), getNumColumns(), mOffset, mValues, getContextPtr() );
    }
    else
    {
        HArray<ValueType> tmpValues;

        DIAUtils::convertDIA2CSR( csrIA, csrJA, tmpValues,
                                  getNumRows(), getNumColumns(), mOffset, mValues, getContextPtr() );

        HArrayUtils::_assign( csrValues, tmpValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setCSRData(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const _HArray& values )
{
    if ( values.getValueType() == getValueType() )
    {
        // call directly setCSRDataImpl, no conversion but cast

        setCSRDataImpl( numRows, numColumns, ia, ja, 
                        static_cast<const HArray<ValueType>&>( values ) );
    }
    else
    {
        // call setCSRDataImpl with converted values

        setCSRDataImpl( numRows, numColumns, ia, ja, 
                        utilskernel::convertHArray<ValueType>( values, getContextPtr() ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values )
{
    SCAI_REGION( "Storage.DIA.setCSR" )

    SCAI_LOG_INFO( logger, "setCSRData " << numRows << " x " << numColumns 
                            << ", ia = " << ia << ", ja = " << ja << ", values = " << values )

    IndexType numValues = ja.size();

    if ( ia.size() == numRows )
    {
        HArray<IndexType> tmpOffsets;
        IndexType total = CSRUtils::sizes2offsets( tmpOffsets, ia, getContextPtr() );
        SCAI_ASSERT_EQUAL( total, numValues, "sizes do not sum up correctly" )
        setCSRDataImpl( numRows, numColumns, tmpOffsets, ja, values );
        return;
    }

    SCAI_ASSERT_EQ_DEBUG( numRows + 1, ia.size(), "illegal CSR ia offset array" )
    SCAI_ASSERT_EQ_DEBUG( ja.size(), values.size(), "serious size mismatch for CSR arrays." )

    SCAI_ASSERT_DEBUG( CSRUtils::validOffsets( ia, numValues, getContextPtr() ), "illegal CSR offset array" );

    SCAI_ASSERT_DEBUG( HArrayUtils::validIndexes( ja, numColumns, getContextPtr() ),
                       "CSR ja array contains illegal column indexes, #columns = " << numColumns );

    _MatrixStorage::setDimension( numRows, numColumns );

    DIAUtils::convertCSR2DIA( mOffset, mValues, numRows, numColumns, ia, ja, values, getContextPtr() );
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

    return DIAUtils::maxNorm( getNumRows(), getNumColumns(), mOffset, mValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType DIAStorage<ValueType>::getValue( const IndexType i, const IndexType j ) const
{
    IndexType pos = DIAUtils::getValuePos( i, j, getNumRows(), getNumColumns(), mOffset, getContextPtr() );

    ValueType val = 0;

    if ( pos != invalidIndex )
    {
        SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mValues.size(), "illegal value position for ( " << i << ", " << j << " )" );

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

    IndexType pos = DIAUtils::getValuePos( i, j, getNumRows(), getNumColumns(), mOffset, getContextPtr() );

    if ( pos == invalidIndex )
    {
        COMMON_THROWEXCEPTION( "DIA storage has no entry ( " << i << ", " << j << " ) " )
    }

    SCAI_ASSERT_VALID_INDEX_DEBUG( pos, mValues.size(), "illegal value position for ( " << i << ", " << j << " )" );

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
size_t DIAStorage<ValueType>::getMemoryUsage() const
{
    size_t memoryUsage = _MatrixStorage::_getMemoryUsage();

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

    return DIAUtils::gemv( result, alpha, x, beta, y,
                           getNumRows(), getNumColumns(), mOffset, mValues,
                           op, async, getContextPtr() );
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
    SCAI_REGION( "Storage.DIA.jacobiIterate" )

    SCAI_LOG_INFO( logger, *this << ": Jacobi iteration for local matrix data." )

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of solution and oldSolution unsupported" )
    }

    // matrix must be square

    SCAI_ASSERT_EQ_DEBUG( getNumRows(), getNumColumns(), "jacobi iteration step only on square matrix storage" )

    bool async = false;  // no sync token, NULL return can be ignored

    DIAUtils::jacobi( solution, omega, oldSolution, rhs, 
                      getNumRows(), mOffset, mValues, async, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAStorage<ValueType>::jacobiIterateAsync(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.DIA.jacobiIterateAsync" )

    bool async = true;  // call will return valid SyncToken

    SyncToken* token = DIAUtils::jacobi( solution, omega, oldSolution, rhs,
                           getNumRows(), mOffset, mValues, async, getContextPtr() );

    if ( token == NULL )
    {
        // there was no asynchronous execution at all

        token = new tasking::NoSyncToken();
    }

    return token;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void DIAStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& solution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldSolution,
    const ValueType omega ) const
{
    SCAI_REGION( "Storage.DIA.jacobiIterateHalo" )

    // matrix must be square

    DIAUtils::jacobiHalo( solution, omega, localDiagonal, oldSolution, 
                          getNumRows(), getNumColumns(), mOffset, mValues, getContextPtr() );
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

} /* end namespace lama */

} /* end namespace scai */
