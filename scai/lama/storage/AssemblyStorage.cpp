/**
 * @file AssemblyStorage.cpp
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
 * @brief AssemblyStorage.cpp
 * @author Jiri Kraus, Thomas Brandes
 * @date 07.11.2011
 */

// hpp
#include <scai/lama/storage/AssemblyStorage.hpp>
#include <scai/lama/storage/CSRStorage.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRUtils.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>
#include <scai/hmemo.hpp>
#include <scai/common/macros/print_string.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>
#include <scai/common/Constants.hpp>
#include <scai/common/macros/instantiate.hpp>

namespace scai
{

using namespace hmemo;

using utilskernel::HArrayUtils;

namespace lama
{

/* --------------------------------------------------------------------------- */

SCAI_LOG_DEF_TEMPLATE_LOGGER( template<typename ValueType>, AssemblyStorage<ValueType>::logger,
                              "MatrixStorage.AssemblyStorage" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>::AssemblyStorage() :

    MatrixStorage<ValueType>( 0, 0, hmemo::Context::getHostPtr() ), 
    mRows( 0 ), 
    mNumValues( 0 )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>::AssemblyStorage( const AssemblyStorage<ValueType>& other ) : 

    MatrixStorage<ValueType>( other.getNumRows(), other.getNumColumns(), other.getContextPtr() ), 
    mRows( other.mRows ), 
    mNumValues( other.mNumValues )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>::AssemblyStorage( const _MatrixStorage& other ) : 

    MatrixStorage<ValueType>( other.getNumRows(), other.getNumColumns(), hmemo::Context::getHostPtr() )

{
    assign( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::assign( const _MatrixStorage& other )
{
    // translate virtual call to specific template call via wrapper

    mepr::StorageWrapper<AssemblyStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::assignImpl( this, other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void AssemblyStorage<ValueType>::assignImpl( const MatrixStorage<OtherValueType>& other )
{
    if ( static_cast<const _MatrixStorage*>( &other ) == this )
    {
        SCAI_LOG_INFO( logger, typeName() << ": self assign, skipped, matrix = " << other )
    }
    else
    {
        HArray<IndexType> csrIA;
        HArray<IndexType> csrJA;
        HArray<ValueType> csrValues;

        other.buildCSRData( csrIA, csrJA, csrValues );

        // setCSRDataImpl can deal with alias and takes advantage of it 
        // and it will set all relevant attributes of this storage correctly.

        setCSRDataImpl( other.getNumRows(), other.getNumColumns(), csrIA, csrJA, csrValues );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>::AssemblyStorage(
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow /* = 10*/ ) : 

    MatrixStorage<ValueType>( numRows, numColumns, hmemo::Context::getHostPtr() ), 
    mRows( numRows ), 
    mNumValues( 0 )
{
    SCAI_LOG_INFO( logger,
                   "Creating with " << getNumRows() << " rows, " << getNumColumns() << " columns, " << numValuesPerRow << " values per row." )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        SCAI_LOG_TRACE( logger, "Reserving storage for row " << i )
        mRows[i].reserve( numValuesPerRow );
    }

    SCAI_LOG_DEBUG( logger, "Created." )
}

template<typename ValueType>
AssemblyStorage<ValueType>::~AssemblyStorage()
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>& AssemblyStorage<ValueType>::operator=(
    const AssemblyStorage<ValueType>& other )
{
    _MatrixStorage::setDimension( other.getNumRows(), other.getNumColumns() );
    mRows = other.mRows;
    mNumValues = other.mNumValues;
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::swap( _MatrixStorage& other )
{
    SCAI_ASSERT_EQ_ERROR( getFormat(), other.getFormat(), "swap only for same storage format" )
    SCAI_ASSERT_EQ_ERROR( this->getValueType(), other.getValueType(), "swap only for same value type" )

    // only in debug mode use the more expensive dynamic cast for verification

    SCAI_ASSERT_DEBUG( dynamic_cast<AssemblyStorage<ValueType>* >( &other ), "illegal storage to swap" )

    swapImpl( reinterpret_cast<AssemblyStorage<ValueType>& >( other ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::swapImpl( AssemblyStorage<ValueType>& other )
{
    std::swap( mNumValues, other.mNumValues );
    mRows.swap( other.mRows );
    _MatrixStorage::swap( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>& AssemblyStorage<ValueType>::operator=( const _MatrixStorage& other )
{
    assign( other );
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::allocate( const IndexType numRows, const IndexType numColumns )
{
    SCAI_LOG_INFO( logger, "allocate sparse assembly storage " << numRows << " x " << numColumns )
    _MatrixStorage::setDimension( numRows, numColumns );
    mNumValues = 0;
    mRows.clear();
    mRows.resize( numRows );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
Format AssemblyStorage<ValueType>::getFormat() const
{
    return Format::ASSEMBLY;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::prefetch( const ContextPtr /* location */ ) const
{
    // not supported, no implementation required
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::wait() const
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::clear()
{
    // clear is here like a purge due to destructor of mRows
    _MatrixStorage::setDimension( 0, 0 );
    mNumValues = 0;
    mRows.clear();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::purge()
{
    _MatrixStorage::setDimension( 0, 0 );
    mNumValues = 0;
    mRows.clear();
}

/* --------------------------------------------------------------------------- */

#ifdef SCAI_ASSERT_LEVEL_OFF
template<typename ValueType>
void AssemblyStorage<ValueType>::check( const char* ) const
{}
#else
template<typename ValueType>
void AssemblyStorage<ValueType>::check( const char* msg ) const
{
    SCAI_ASSERT_EQ_ERROR( getNumRows(), static_cast<IndexType>( mRows.size() ), msg << ", serious mismatch" )

    // TODO check that numValues is equal sum( mRows[i].ja.size() ), 0 <= i < getNumRows()
    // TODO check for valid column indexes
}
#endif

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>* AssemblyStorage<ValueType>::newMatrixStorage( const IndexType numRows, const IndexType numColumns ) const
{
    return new AssemblyStorage<ValueType>( numRows, numColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> AssemblyStorage<ValueType>::l1Norm() const
{
    RealType<ValueType> val = 0;

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        for ( size_t jj = 0; jj < mRows[i].values.size(); ++jj )
        {
            val += common::Math::abs( mRows[i].values[jj] );
        }
    }

    return val;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> AssemblyStorage<ValueType>::l2Norm() const
{
    RealType<ValueType> val = 0;

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        for ( size_t jj = 0; jj < mRows[i].values.size(); ++jj )
        {
            auto tmp = common::Math::abs( mRows[i].values[jj] );
            val += tmp * tmp;
        }
    }

    return common::Math::sqrt( val );
}


/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> AssemblyStorage<ValueType>::maxNorm() const
{
    RealType<ValueType> maxval = 0;

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        const std::vector<ValueType>& values = mRows[i].values;

        for ( size_t jj = 0; jj < values.size(); ++jj )
        {
            auto val = common::Math::abs( mRows[i].values[jj] );

            maxval = common::Math::max( common::Math::abs( val ), maxval );
        }
    }

    return maxval;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
size_t AssemblyStorage<ValueType>::getMemoryUsageImpl() const
{
    size_t memoryUsage = 0;
    memoryUsage += sizeof( IndexType ) * mNumValues;
    memoryUsage += sizeof( ValueType ) * mNumValues;
    return memoryUsage;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
AssemblyStorage<ValueType>::Row::Row()
{
}

template<typename ValueType>
AssemblyStorage<ValueType>::Row::Row( const IndexType numValuesPerRow )
{
    ja.reserve( numValuesPerRow );
    values.reserve( numValuesPerRow );
}

template<typename ValueType>
AssemblyStorage<ValueType>::Row::Row( const typename AssemblyStorage<ValueType>::Row& other )
    : ja( other.ja ), values( other.values )
{
}

template<typename ValueType>
typename AssemblyStorage<ValueType>::Row& AssemblyStorage<ValueType>::Row::operator=(
    const typename AssemblyStorage<ValueType>::Row& other )
{
    ja = other.ja;
    values = other.values;
    return *this;
}

template<typename ValueType>
void AssemblyStorage<ValueType>::Row::reserve( const IndexType numValuesPerRow )
{
    ja.reserve( numValuesPerRow );
    values.reserve( numValuesPerRow );
}

template<typename ValueType>
void AssemblyStorage<ValueType>::Row::scale( const ValueType val )
{
    for ( size_t i = 0; i < values.size(); i++ )
    {
        values[i] *= val;
    }
}

template<typename ValueType>
void AssemblyStorage<ValueType>::Row::conj()
{
    for ( size_t i = 0; i < values.size(); i++ )
    {
        values[i] = common::Math::conj( values[i] );
    }
}

template<typename ValueType>
ValueType AssemblyStorage<ValueType>::operator()( const IndexType i, const IndexType j ) const
{
    if ( j >= getNumColumns() )
    {
        COMMON_THROWEXCEPTION( "Passed column Index " << j << " exceeds column count " << getNumColumns() << "." )
    }

    const std::vector<IndexType>& rJA = mRows[i].ja;

    for ( size_t k = 0; k < mRows[i].ja.size(); ++k )
    {
        if ( j == rJA[k] )
        {
            const std::vector<ValueType>& rValues = mRows[i].values;
            return rValues[k];
        }
    }

    return static_cast<ValueType>( 0.0 );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::print( std::ostream& stream ) const
{
    using std::endl;
    stream << "AssemblyStorage " << getNumRows() << " x " << getNumColumns() << ", #values = " << mNumValues << endl;

    for ( IndexType i = 0; i < getNumRows(); i++ )
    {
        const Row& row = mRows[i];
        stream << "Row " << i << " ( " << row.ja.size() << " ) :";

        for ( size_t k = 0; k < row.ja.size(); ++k )
        {
            stream << " " << row.ja[k] << ":" << row.values[k];
        }

        stream << endl;
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setValue(
    const IndexType i,
    const IndexType j,
    const ValueType val,
    const common::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX( i, getNumRows(), "illegal row index" )
    SCAI_ASSERT_VALID_INDEX( j, getNumColumns(), "illegal col index" )

    {
        const std::vector<IndexType>& rJA = mRows[i].ja;

        for ( size_t k = 0; k < mRows[i].ja.size(); ++k )
        {
            if ( j == rJA[k] )
            {
                std::vector<ValueType>& wValues = mRows[i].values;
                ValueType* opnd = &wValues[k];
                utilskernel::OpenMPUtils::setVal( opnd, 1, val, op );
                return;
            }
        }
    }

    ValueType newValue = 0;

    utilskernel::OpenMPUtils::setVal( &newValue, 1, val, op );

    SCAI_LOG_TRACE( logger, "set( " << i << ", " << j << ", " << val << ") : new entry " )
    std::vector<IndexType>& wJA = mRows[i].ja;
    std::vector<ValueType>& wValues = mRows[i].values;
    wValues.push_back( newValue );
    wJA.push_back( j );
    ++mNumValues;

    // if we have a diagonal element and the row was not empty before we need to
    // fix the diagonal property

    if ( i == j && ( wJA.size() - 1 > 0 ) )
    {
        SCAI_LOG_TRACE( logger, "diagonal element swapped to first element of row" )
        std::swap( wValues[0], wValues[wValues.size() - 1] );
        std::swap( wJA[0], wJA[wJA.size() - 1] );
    }
}

template<typename ValueType>
void AssemblyStorage<ValueType>::set( const IndexType i, const IndexType j, const ValueType value )
{
    setValue( i, j, value );
}

template<typename ValueType>
const std::vector<IndexType>&
AssemblyStorage<ValueType>::getJa( const IndexType i ) const
{
    return mRows[i].ja;
}

template<typename ValueType>
std::vector<IndexType>& AssemblyStorage<ValueType>::getJa( const IndexType i )
{
    return mRows[i].ja;
}

template<typename ValueType>
const std::vector<ValueType>& AssemblyStorage<ValueType>::getValues( const IndexType i ) const
{
    return mRows[i].values;
}

template<typename ValueType>
IndexType AssemblyStorage<ValueType>::getNumValues() const
{
    return mNumValues;
}

template<typename ValueType>
void AssemblyStorage<ValueType>::setSparseRow(
    const IndexType i,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values )
{
    //SCAI_ASSERT_EQUAL_ERROR( ja.size(), values.size() )
    #pragma omp atomic
    mNumValues -= mRows[i].ja.size();
    mRows[i].ja.resize( ja.size() );
    mRows[i].values.resize( values.size() );
    ReadAccess<IndexType> rJA( ja );
    ReadAccess<ValueType> rValues( values );

    for ( IndexType k = 0; k < ja.size(); ++k )
    {
        mRows[i].ja[k] = rJA[k];
        mRows[i].values[k] = rValues[k];
    }

    const IndexType rowSize = mRows[i].ja.size();

    #pragma omp atomic
    mNumValues += rowSize;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::fixDiagonalProperty( const IndexType i )

{
    // fix diagonal property if necessary
    if ( i >= getNumColumns() )
    {
        return;
    }

    if ( mRows[i].ja.size() == 0 )
    {
        #pragma omp atomic
        ++mNumValues;
        mRows[i].ja.push_back( i );
        mRows[i].values.push_back( static_cast<ValueType>( 0.0 ) );
        return;
    }

    if ( mRows[i].ja[0] == i )
    {
        return;
    }

    // try to find diagonal element
    std::vector<IndexType>& wJA = mRows[i].ja;
    std::vector<ValueType>& wValues = mRows[i].values;

    for ( size_t k = 0; k < wJA.size(); ++k )
    {
        if ( i == wJA[k] )
        {
            std::swap( wValues[0], wValues[k] );
            std::swap( wJA[0], wJA[k] );
            break;
        }
    }

    if ( wJA[0] != i )
    {
        #pragma omp atomic
        ++mNumValues;
        wJA.push_back( i );
        wValues.push_back( static_cast<ValueType>( 0.0 ) );
        std::swap( wValues[0], wValues[wValues.size() - 1] );
        std::swap( wJA[0], wJA[wValues.size() - 1] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setNumColumns( const IndexType numColumns )
{
    _MatrixStorage::setDimension( getNumRows(), numColumns );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setIdentity( const IndexType n )
{
    allocate( n, n );

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        set( i, i, static_cast<ValueType>( 1.0 ) );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void AssemblyStorage<ValueType>::setCSRDataImpl(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<OtherValueType>& values )
{
    ContextPtr loc = getContextPtr();

    if ( ia.size() == numRows )
    {
        // offset array required
        HArray<IndexType> offsets;
        sparsekernel::CSRUtils::sizes2offsets( offsets, ia, loc );
        setCSRDataImpl( numRows, numColumns, offsets, ja, values );
        return;
    }

    IndexType numValues = ja.size();

    SCAI_ASSERT_EQUAL( values.size(), numValues, "size misamtch" );
    SCAI_ASSERT_EQUAL( ia.size(), numRows + 1, "size misamtch" );
    SCAI_ASSERT( HArrayUtils::isSorted( ia, common::CompareOp::LE, loc ),
                 "illegal offset array, not ascending entries" );
    SCAI_ASSERT_EQUAL( HArrayUtils::getVal( ia, numRows ), numValues,
                       "illegal offset array, not ascending entries" );
    SCAI_ASSERT( HArrayUtils::validIndexes( ja, numColumns ),
                 "illegal column indexes, #colums = " << numColumns );
    // no more error checks here on the sizes, but on the content
    ReadAccess<IndexType> csrIA( ia );
    ReadAccess<IndexType> csrJA( ja );
    ReadAccess<OtherValueType> csrValues( values );
    _MatrixStorage::setDimension( numRows, numColumns );
    mNumValues = numValues;
    mRows.resize( getNumRows() );
    SCAI_ASSERT_EQUAL_ERROR( csrIA[numRows], numValues )
    SCAI_LOG_DEBUG( logger, "fill " << *this << " with csr data, " << numValues << " non-zero values" )

    for ( IndexType i = 0; i < numRows; ++i )
    {
        const IndexType n = csrIA[i + 1] - csrIA[i];
        Row& row = mRows[i];
        row.ja.resize( n );
        row.values.resize( n );
        IndexType offset = 0;

        for ( IndexType k = csrIA[i]; k < csrIA[i + 1]; ++k )
        {
            row.ja[offset] = csrJA[k];
            row.values[offset] = static_cast<ValueType>( csrValues[k] );
            ++offset;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void AssemblyStorage<ValueType>::setDIADataImpl(
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
void AssemblyStorage<ValueType>::buildCSRSizes( hmemo::HArray<IndexType>& ia ) const
{
    auto csrIA = hostWriteOnlyAccess( ia, getNumRows() + 1 );

    // build csrSizes in ia

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        csrIA[i] = static_cast<IndexType>( mRows[i].ja.size() );
    }

    csrIA.resize( getNumRows() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void AssemblyStorage<ValueType>::buildCSR(
    HArray<IndexType>& ia,
    HArray<IndexType>* ja,
    HArray<OtherValueType>* values,
    const ContextPtr /* loc */ ) const
{
    SCAI_LOG_INFO( logger, *this << ": build CSR data from it" )

    {
        auto csrIA = hostWriteOnlyAccess( ia, getNumRows() + 1 );

        // build csrSizes in ia

        for ( IndexType i = 0; i < getNumRows(); ++i )
        {
            csrIA[i] = static_cast<IndexType>( mRows[i].ja.size() );
        }
    
        csrIA.resize( getNumRows() );
    }

    if ( ja == NULL || values == NULL )
    {
        return;
    }

    HArrayUtils::scan1( ia );

    auto csrIA = hostReadAccess( ia );
    auto csrJA = hostWriteOnlyAccess( *ja, mNumValues );
    auto csrValues = hostWriteOnlyAccess( *values, mNumValues );

    #pragma omp parallel for
    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        IndexType offset = 0;

        for ( IndexType k = csrIA[i]; k < csrIA[i + 1]; ++k )
        {
            csrJA[k] = mRows[i].ja[offset];
            csrValues[k] = static_cast<OtherValueType>( mRows[i].values[offset] );
            ++offset;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setCSRData(
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const hmemo::_HArray& values )
{
    mepr::StorageWrapper<AssemblyStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
        setCSRDataImpl( this, numRows, numColumns, ia, ja, values );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::buildCSRData(
    hmemo::HArray<IndexType>& csrIA,
    hmemo::HArray<IndexType>& csrJA,
    hmemo::_HArray& csrValues ) const
{
    mepr::StorageWrapper<AssemblyStorage, SCAI_NUMERIC_TYPES_HOST_LIST>::
        buildCSRDataImpl( this, csrIA, csrJA, csrValues, getContextPtr() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setDiagonalV( const HArray<ValueType>& diagonal )
{
    IndexType numDiagonalElements = diagonal.size();
    ReadAccess<ValueType> rDiagonal( diagonal );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        set( i, i, rDiagonal[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::getDiagonal( HArray<ValueType>& diagonal ) const
{
    const IndexType numDiagonalElements = std::min( getNumColumns(), getNumRows() );
    WriteOnlyAccess<ValueType> wDiagonal( diagonal, numDiagonalElements );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        wDiagonal[i] = operator()( i, i );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::assignDiagonal( const HArray<ValueType>& diagonal )
{
    const IndexType size = diagonal.size();

    _MatrixStorage::setDimension( size, size );

    mRows.clear();
    mRows.resize( size );
    mNumValues = 0;

    setDiagonalV( diagonal );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::getRow( HArray<ValueType>& row, const IndexType i ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    auto wRow = hostWriteOnlyAccess( row, getNumColumns() );

    for ( IndexType j = 0; j < getNumColumns(); ++j )
    {
        wRow[j] = 0;
    }

    const std::vector<IndexType>& ja = mRows[i].ja;
    const std::vector<ValueType>& values = mRows[i].values;

    SCAI_ASSERT_EQUAL_DEBUG( ja.size(), values.size() )

    for ( size_t k = 0; k < ja.size(); ++k )
    {
        const IndexType j = ja[k];

        SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "illegal col index at row " << i << ", ja[" << k << "]" )
        wRow[j] = values[k];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::getSparseRow( hmemo::HArray<IndexType>& jA, hmemo::HArray<ValueType>& values, const IndexType i ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )

    SCAI_ASSERT_EQ_ERROR( mRows[i].ja.size(), mRows[i].values.size(), "serious mismatch in row " << i )

    jA = HArray<IndexType>( mRows[i].ja );
    values = HArray<ValueType>( mRows[i].values );
}

template<typename ValueType>
void AssemblyStorage<ValueType>::getSparseColumn( hmemo::HArray<IndexType>& iA, hmemo::HArray<ValueType>& values, const IndexType j ) const
{
    IndexType colSize = 0;

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        for ( size_t jj = 0; jj < mRows[i].ja.size(); ++jj )
        {
            if ( mRows[i].ja[jj] == j )
            {
                colSize++;
                break;
            }
        }
    }

    auto wIA = hostWriteOnlyAccess( iA, colSize );
    auto wValues = hostWriteOnlyAccess( values, colSize );

    colSize = 0;  // reset counter for column entries

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        for ( size_t jj = 0; jj < mRows[i].ja.size(); ++jj )
        {
            if ( mRows[i].ja[jj] == j )
            {
                wIA[colSize] = i;
                wValues[colSize] = mRows[i].values[jj];
                colSize++;
                break;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::getColumn( HArray<ValueType>& column, const IndexType j ) const
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( j, getNumColumns(), "column index out of range" )

    auto wColumn = hostWriteOnlyAccess( column, getNumRows() );

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        wColumn[i] = 0;

        for ( size_t jj = 0; jj < mRows[i].ja.size(); ++jj )
        {
            if ( mRows[i].ja[jj] == j )
            {
                wColumn[i] = mRows[i].values[jj];
                break;
            }
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setRow( const HArray<ValueType>& row, const IndexType i, const common::BinaryOp op )
{
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, getNumRows(), "row index out of range" )
    SCAI_ASSERT_GE_DEBUG( row.size(), getNumColumns(), "row array to small for set" )

    // ToDo write more efficient kernel routine for setting a row

    auto rRow = hostReadAccess( row );

    for ( IndexType j = 0; j < getNumColumns(); ++j )
    {
        if ( rRow[j] == common::Constants::ZERO )
        {
            continue;
        }

        setValue( i, j, rRow[j], op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setColumn( const HArray<ValueType>& column, const IndexType j, const common::BinaryOp op )
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

        setValue( i, j, rColumn[i], op );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::setDiagonal( const ValueType value )
{
    const IndexType numDiagonalElements = common::Math::min( getNumColumns(), getNumRows() );

    for ( IndexType i = 0; i < numDiagonalElements; ++i )
    {
        set( i, i, value );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::scale( const ValueType value )
{
    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        mRows[i].scale( value );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::conj()
{
    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        mRows[i].conj();
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::scaleRows( const HArray<ValueType>& diagonal )
{
    SCAI_ASSERT_EQ_ERROR( getNumRows(), diagonal.size(), "not one element for each row" )

    ReadAccess<ValueType> rDiagonal( diagonal );

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        mRows[i].scale( rDiagonal[i] );
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::scaleColumns( const HArray<ValueType>& diagonal )
{
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), diagonal.size(), "not one element for each col" )

    ReadAccess<ValueType> rDiagonal( diagonal );

    for ( IndexType i = 0; i < getNumRows(); ++i )
    {
        for ( size_t jj = 0 ; jj < mRows[i].ja.size(); ++jj )
        {
            mRows[i].values[jj] *= rDiagonal[mRows[i].ja[jj]];
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::matrixTimesVector(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    static bool warn = true;   

    if ( warn )
    {
        SCAI_LOG_WARN( logger, "matrixTimesVector<" << getValueType() 
                                << "> not supported on an assembly storage, will be converted to CSR" )

        warn = false;  // print warning only once
    }

    CSRStorage<ValueType> csr;
    csr.assign( *this );
    csr.matrixTimesVector( result, alpha, x, beta, y, op );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void AssemblyStorage<ValueType>::writeAt( std::ostream& stream ) const
{
    stream << "AssemblyStorage<" << common::getScalarType<ValueType>() << ">("
           << " size = " << getNumRows() << " x " << getNumColumns()
           << ", #values = " << mNumValues << " )";
}

template<typename ValueType>
_MatrixStorage* AssemblyStorage<ValueType>::create()
{
    return new AssemblyStorage<ValueType>();
}

template<typename ValueType>
MatrixStorageCreateKeyType AssemblyStorage<ValueType>::createValue()
{
    return MatrixStorageCreateKeyType( Format::ASSEMBLY, common::getScalarType<ValueType>() );
}

template<typename ValueType>
std::string AssemblyStorage<ValueType>::initTypeName()
{
    std::stringstream s;
    s << std::string( "AssemblyStorage<" ) << common::getScalarType<ValueType>() << std::string( ">" );
    return s.str();
}

template<typename ValueType>
const char* AssemblyStorage<ValueType>::typeName()
{
    static const std::string s = initTypeName();
    return  s.c_str();
}

template<typename ValueType>
const char* AssemblyStorage<ValueType>::getTypeName() const
{
    return typeName();
}

template<typename ValueType>
MatrixStorageCreateKeyType AssemblyStorage<ValueType>::getCreateValue() const
{
    return MatrixStorageCreateKeyType( Format::ASSEMBLY, common::getScalarType<ValueType>() );
}

/* ========================================================================= */
/*       Template specializattions and instantiations                        */
/* ========================================================================= */

SCAI_COMMON_INST_CLASS( AssemblyStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
