/**
 * @file MatrixStorage.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of methods for common base class of all matrix storage formats.
 * @author Thomas Brandes
 * @date 27.04.2011
 * @since 1.0.0
 */

// hpp
#include <scai/lama/storage/MatrixStorage.hpp>

// others
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/StorageMethods.hpp>

#include <scai/lama/distribution/Distribution.hpp>
#include <scai/lama/distribution/Redistributor.hpp>
#include <scai/lama/distribution/Halo.hpp>

#include <scai/lama/LAMAInterface.hpp>
#include <scai/lama/StorageIO.hpp>

#include <scai/lama/openmp/OpenMPUtils.hpp>
#include <scai/lama/openmp/OpenMPCSRUtils.hpp>

#include <scai/lama/LAMATypes.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

// tracing
#include <scai/tracing.hpp>


#include <scai/lama/LAMATypes.hpp>

// boost
#include <scai/common/bind.hpp>
#include <boost/preprocessor.hpp>
#include <boost/assign/list_of.hpp>
#include <boost/unordered_map.hpp>

using scai::tasking::SyncToken;
using scai::tasking::TaskSyncToken;

using namespace scai::hmemo;

namespace scai
{

namespace lama
{

SCAI_LOG_DEF_LOGGER( _MatrixStorage::logger, "MatrixStorage" )

_MatrixStorage::_MatrixStorage()

    : mNumRows( 0 ), mNumColumns( 0 ), mRowIndexes(), mCompressThreshold( 0.0f ), mDiagonalProperty(
        false ), mContext( Context::getContextPtr( context::Host ) )
{
    SCAI_LOG_DEBUG( logger, "constructed MatrixStorage()" )
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::setDimension( const IndexType numRows, const IndexType numColumns )
{
    // in any case set dimensions
    mNumRows = numRows;
    mNumColumns = numColumns;
    // due to new settings assume that diagonalProperty, rowIndexes become invalid
    mDiagonalProperty = false;
    mRowIndexes.clear();
    // but do not reset threshold
}

/* ---------------------------------------------------------------------------------- */

_MatrixStorage::~_MatrixStorage()
{
    SCAI_LOG_DEBUG( logger, "_MatrixStorage" )
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::writeAt( std::ostream& stream ) const
{
    stream << " MatrixStorage: (" << mNumRows << " x " << mNumColumns << ")";
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::setCompressThreshold( float ratio )
{
    if ( ratio < 0.0 || ratio > 1.0 )
    {
        COMMON_THROWEXCEPTION( "Illegal threshold " << ratio << ", must be from 0.0 to 1.0" )
    }

    mCompressThreshold = ratio;
    SCAI_LOG_INFO( logger, "set compress threshold, ratio = " << ratio << " : " << *this )
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::swap( _MatrixStorage& other )
{
    std::swap( mNumRows, other.mNumRows );
    std::swap( mNumColumns, other.mNumColumns );
    mRowIndexes.swap( other.mRowIndexes );
    std::swap( mCompressThreshold, other.mCompressThreshold );
    std::swap( mDiagonalProperty, other.mDiagonalProperty );
    std::swap( mContext, other.mContext );
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::_assignTranspose( const _MatrixStorage& other )
{
    mNumRows = other.mNumColumns;
    mNumColumns = other.mNumRows;
    mRowIndexes.clear();
    mCompressThreshold = other.mCompressThreshold;
    mDiagonalProperty = false;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::_assign( const _MatrixStorage& other )
{
    mNumRows = other.mNumRows;
    mNumColumns = other.mNumColumns;
    mRowIndexes.clear();
    mCompressThreshold = other.mCompressThreshold;
    mDiagonalProperty = false;
}

/* --------------------------------------------------------------------------- */

_MatrixStorage& _MatrixStorage::operator=( const _MatrixStorage& other )
{
    assign( other ); // assign can deal with all kind of storage formats/types
    return *this;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::resetDiagonalProperty()
{
    mDiagonalProperty = checkDiagonalProperty();
    SCAI_LOG_DEBUG( logger, *this << ": diagonal property = " << mDiagonalProperty )
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::setContext( ContextPtr context )
{
    if ( context.get() != mContext.get() )
    {
        SCAI_LOG_DEBUG( logger, *this << ": new location = " << *context << ", old location = " << *mContext )
    }

    mContext = context;
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::localize( const _MatrixStorage& global, const Distribution& rowDist )
{
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), global.getNumColumns() )
    SCAI_ASSERT_EQUAL_ERROR( global.getNumRows(), rowDist.getGlobalSize() )
    COMMON_THROWEXCEPTION( "No default implementation for localize available, matrix = " << *this )
}

/* --------------------------------------------------------------------------- */

IndexType _MatrixStorage::getNumValues() const
{
    // Default implementation builds sum of row sizes
    LAMAArray<IndexType> sizes;
    buildCSRSizes( sizes );
    ReadAccess<IndexType> csrSizes( sizes );
    IndexType numValues = OpenMPUtils::sum( csrSizes.get(), mNumRows );
    return numValues;
}

/* ---------------------------------------------------------------------------------- */

static const boost::unordered_map<MatrixStorageFormat, const char*> formatToString = boost::assign::map_list_of(
            Format::CSR, "CSR" )( Format::ELL, "ELL" )( Format::DIA, "DIA" )( Format::JDS, "JDS" )( Format::COO,
                    "COO" )(
            Format::DENSE, "DENSE" )( Format::ASSEMBLY, "ASSEMBLY" )( Format::UNDEFINED, "UNDEFINED" );

const char* format2Str( const MatrixStorageFormat storageFormat )
{
    return formatToString.at( storageFormat );
}

MatrixStorageFormat str2Format( const char* str )
{
    for ( int format = Format::CSR; format < Format::UNDEFINED; ++format )
    {
        if ( strcmp( formatToString.at( MatrixStorageFormat( format ) ), str ) == 0 )
        {
            return MatrixStorageFormat( format );
        }
    }

    return Format::UNDEFINED;
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::offsets2sizes( LAMAArray<IndexType>& offsets )
{
    const IndexType n = offsets.size() - 1;
    WriteAccess<IndexType> writeSizes( offsets );

    // the following loop  is not parallel

    for ( IndexType i = 0; i < n; i++ )
    {
        writeSizes[i] = writeSizes[i + 1] - writeSizes[i];
    }

    writeSizes.resize( n );
}

/* ---------------------------------------------------------------------------------- */

void _MatrixStorage::offsets2sizes( LAMAArray<IndexType>& sizes, const LAMAArray<IndexType>& offsets )
{
    if ( &sizes == &offsets )
    {
        SCAI_LOG_WARN( logger, "offset2sizes: sizes and offsets are same array" )
        offsets2sizes( sizes );
        return;
    }

    const IndexType n = offsets.size() - 1;

    ReadAccess<IndexType> readOffsets( offsets );

    WriteAccess<IndexType> writeSizes( sizes );

    writeSizes.clear(); // old values are not used

    writeSizes.resize( n );

    for ( IndexType i = 0; i < n; i++ )
    {
        writeSizes[i] = readOffsets[i + 1] - readOffsets[i];
    }
}

/* ---------------------------------------------------------------------------------- */

IndexType _MatrixStorage::sizes2offsets( LAMAArray<IndexType>& sizes )
{
    IndexType n = sizes.size();
    WriteAccess<IndexType> writeOffsets( sizes );
    writeOffsets.resize( n + 1 );
    return OpenMPCSRUtils::sizes2offsets( writeOffsets, n );
}

/* ---------------------------------------------------------------------------------- */

size_t _MatrixStorage::getMemoryUsage() const
{
    size_t memoryUsage = 0;
    memoryUsage += 2 * sizeof( IndexType );
    memoryUsage += sizeof( bool );
    memoryUsage += sizeof( float );
    memoryUsage += sizeof( IndexType ) * mRowIndexes.size();
    memoryUsage += getMemoryUsageImpl();
    SCAI_LOG_DEBUG( logger, *this << ": used memory = " << memoryUsage )
    return memoryUsage;
}

template<typename ValueType>
MatrixStorage<ValueType>::MatrixStorage()

    : _MatrixStorage(), mEpsilon( 0 )
{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorage<ValueType>::MatrixStorage( const IndexType numRows, const IndexType numColumns )

    : _MatrixStorage(), mEpsilon( 0 )
{
    setDimension( numRows, numColumns );
    SCAI_LOG_DEBUG( logger, "constructed MatrixStorage<ValueType> for " << numRows << " x " << numColumns << " matrix" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorage<ValueType>::~MatrixStorage()
{
    SCAI_LOG_DEBUG( logger, "~MatrixStorage" )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType MatrixStorage<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::swap( MatrixStorage<ValueType>& other )
{
    _MatrixStorage::swap( other );
    std::swap( mEpsilon, other.mEpsilon );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::convertCSR2CSC(
    LAMAArray<IndexType>& colIA,
    LAMAArray<IndexType>& colJA,
    LAMAArray<ValueType>& colValues,
    const IndexType numColumns,
    const LAMAArray<IndexType>& rowIA,
    const LAMAArray<IndexType>& rowJA,
    const LAMAArray<ValueType>& rowValues,
    const ContextPtr loc )
{
    // ContextPtr loc = Context::getContextPtr( context::Host );
    const IndexType numRows = rowIA.size() - 1;
    const IndexType numValues = rowJA.size();
    SCAI_ASSERT_EQUAL_DEBUG( rowJA.size(), rowValues.size() )
    LAMA_INTERFACE_FN_T( convertCSR2CSC, loc, CSRUtils, Transpose, ValueType )
    SCAI_LOG_INFO( logger,
                   "MatrixStorage::CSR2CSC of matrix " << numRows << " x " << numColumns << ", #nnz = " << numValues << " on " << *loc )
    SCAI_REGION( "Storage.CSR2CSC" )
    WriteOnlyAccess<IndexType> cIA( colIA, loc, numColumns + 1 );
    WriteOnlyAccess<IndexType> cJA( colJA, loc, numValues );
    WriteOnlyAccess<ValueType> cValues( colValues, loc, numValues );
    ReadAccess<IndexType> rIA( rowIA, loc );
    ReadAccess<IndexType> rJA( rowJA, loc );
    ReadAccess<ValueType> rValues( rowValues, loc );
    SCAI_CONTEXT_ACCESS( loc )
    convertCSR2CSC( cIA.get(), cJA.get(), cValues.get(), rIA.get(), rJA.get(), rValues.get(), numRows, numColumns,
                    numValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::buildCSCData(
    LAMAArray<IndexType>& colIA,
    LAMAArray<IndexType>& colJA,
    LAMAArray<ValueType>& colValues ) const
{
    LAMAArray<IndexType> rowIA;
    LAMAArray<IndexType> rowJA;
    LAMAArray<ValueType> rowValues;
    buildCSRData( rowIA, rowJA, rowValues );
    ContextPtr loc = Context::getContextPtr( context::Host );
    convertCSR2CSC( colIA, colJA, colValues, mNumColumns, rowIA, rowJA, rowValues, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::assign( const _MatrixStorage& other )
{
    SCAI_REGION( "Storage.assign" )

    if ( &other == this )
    {
        // self assignments might occur during redistributions
        SCAI_LOG_INFO( logger, *this << ": self assign (skipped)" )
        return;
    }

    _MatrixStorage::_assign( other );
    SCAI_LOG_INFO( logger, *this << ": assign ( " << other << " )" )

    if ( other.getFormat() == Format::CSR )
    {
        // CSR storage has more efficient solution: just set CSR data
        other.copyTo( *this );
        return;
    }

    // If the size of other value type is smaller that this value type, it might be better
    // to use the other value type.

    if ( other.getValueType() == common::scalar::FLOAT && getValueType() == common::scalar::DOUBLE )
    {
        other.copyTo( *this );
        return;
    }

    SCAI_LOG_INFO( logger, *this << ": (default) assign ( " << other << " )" )
    LAMAArray<IndexType> csrIA;
    LAMAArray<IndexType> csrJA;
    LAMAArray<ValueType> csrValues;
    other.buildCSRData( csrIA, csrJA, csrValues );
    IndexType numRows = other.getNumRows();
    IndexType numColumns = other.getNumColumns();
    IndexType numValues = csrJA.size();
    SCAI_LOG_DEBUG( logger, "build CSR data " << numRows << " x " << numColumns << ", #nnz = " << numValues )
    setCSRData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
    SCAI_LOG_INFO( logger, "now assigned: " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::copyTo( _MatrixStorage& other ) const
{
    SCAI_REGION( "Storage.copyTo" )
    // If the size of other value type is smaller that this value type, it might be better
    // to use the other value type.
    SCAI_LOG_INFO( logger, *this << ": (default) copyTo ( " << other << " )" )
    LAMAArray<IndexType> csrIA;
    LAMAArray<IndexType> csrJA;
    LAMAArray<ValueType> csrValues;
    buildCSRData( csrIA, csrJA, csrValues );
    SCAI_ASSERT_EQUAL_DEBUG( csrIA.size(), mNumRows + 1 )
    IndexType numValues = csrJA.size();
    SCAI_ASSERT_EQUAL_DEBUG( csrValues.size(), numValues )
    SCAI_LOG_DEBUG( logger, "build CSR data " << mNumRows << " x " << mNumColumns << ", #nnz = " << numValues )
    other.setCSRData( mNumRows, mNumColumns, numValues, csrIA, csrJA, csrValues );
    SCAI_LOG_INFO( logger, "now assigned: " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::assignTranspose( const MatrixStorage<ValueType>& other )
{
    SCAI_REGION( "Storage.assignTranspose" )
    SCAI_LOG_INFO( logger, *this << ": assignTranspose " << other )
    LAMAArray<IndexType> cscIA;
    LAMAArray<IndexType> cscJA;
    LAMAArray<ValueType> cscValues;
    other.buildCSCData( cscIA, cscJA, cscValues );
    // Compressed sparse column data can be used directly to generate the transposed matrix
    // by interpretation as CSR data.
    setCSRData( other.getNumColumns(), other.getNumRows(), cscJA.size(), cscIA, cscJA, cscValues );
    SCAI_LOG_INFO( logger, "now assigned tranposed: " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::joinRows(
    LAMAArray<IndexType>& outSizes,
    LAMAArray<IndexType>& outJA,
    LAMAArray<ValueType>& outValues,
    const IndexType numLocalRows,
    const LAMAArray<IndexType>& rowIndexes,
    const LAMAArray<IndexType>& inSizes,
    const LAMAArray<IndexType>& inJA,
    const LAMAArray<ValueType>& inValues )
{
    SCAI_REGION( "Storage.joinRows" )
    SCAI_LOG_INFO( logger, "join " << numLocalRows << " rows " )
    {
        WriteOnlyAccess<IndexType> sizes( outSizes, numLocalRows );
        ReadAccess<IndexType> rowSizes( inSizes );
        ReadAccess<IndexType> indexes( rowIndexes );

        // initialize counters (Attention: sizes.size() != rowSizes.size())

        for ( int i = 0; i < sizes.size(); ++i )
        {
            sizes[i] = 0;
        }

        // count elements for each row
        for ( int i = 0; i < rowSizes.size(); ++i )
        {
            sizes[indexes[i]] += rowSizes[i];
        }
    }
    // generate offset array for insertion
    LAMAArray<IndexType> IA;
    {
        WriteOnlyAccess<IndexType> offsets( IA, numLocalRows + 1 );
        ReadAccess<IndexType> sizes( outSizes );
        OpenMPUtils::set( offsets.get(), sizes.get(), numLocalRows );
        OpenMPCSRUtils::sizes2offsets( offsets.get(), numLocalRows );
    }
    WriteAccess<IndexType> tmpIA( IA );
    WriteAccess<IndexType> ja( outJA );
    WriteAccess<ValueType> values( outValues );
    ReadAccess<IndexType> rowSizes( inSizes );
    ReadAccess<IndexType> rowJA( inJA );
    ReadAccess<ValueType> rowValues( inValues );
    ReadAccess<IndexType> indexes( rowIndexes );
    // resize data arrays
    ja.resize( rowJA.size() );
    values.resize( rowValues.size() );
    // insert rows
    IndexType dataIndex = 0;

    for ( int i = 0; i < rowSizes.size(); ++i )
    {
        IndexType currentRow = indexes[i];

        for ( int ii = 0; ii < rowSizes[i]; ++ii )
        {
            // insert data at old position 'dataIndex' in row 'currentRow'
            // at new position 'tmpIA[currentRow]'
            ja[tmpIA[currentRow]] = rowJA[dataIndex];
            values[tmpIA[currentRow]] = rowValues[dataIndex];
            dataIndex++;
            tmpIA[currentRow]++;
        }
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorage<ValueType>& MatrixStorage<ValueType>::operator=( const _MatrixStorage& other )
{
    assign( other ); // assign can deal with all kind of storage formats/types
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::joinHalo(
    const _MatrixStorage& localData,
    const _MatrixStorage& haloData,
    const Halo& halo,
    const Distribution& colDist,
    const bool attemptDiagonalProperty )
{
    SCAI_REGION( "Storage.joinHalo" )
    SCAI_LOG_INFO( logger,
                   "join local = " << localData << " with diag = " << localData.hasDiagonalProperty() << " and halo = " << haloData << ", col dist = " << colDist )
    //  Default solution joins storage data via the CSR format
    //  Note: this solution works also for *this == localData or haloData
    LAMAArray<IndexType> localIA;
    LAMAArray<IndexType> localJA;
    LAMAArray<ValueType> localValues;
    localData.buildCSRData( localIA, localJA, localValues );
    SCAI_LOG_DEBUG( logger, "local CSR: ia = " << localIA << ", ja = " << localJA << ", values = " << localValues )
    // map back the local indexes to global column indexes
    {
        IndexType numValues = localJA.size();
        WriteAccess<IndexType> ja( localJA );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            ja[i] = colDist.local2global( ja[i] );
        }
    }
    LAMAArray<IndexType> haloIA;
    LAMAArray<IndexType> haloJA;
    LAMAArray<ValueType> haloValues;
    haloData.buildCSRData( haloIA, haloJA, haloValues );
    SCAI_LOG_DEBUG( logger, "halo CSR: ia = " << haloIA << ", ja = " << haloJA << ", values = " << haloValues )
    // map back the halo indexes to global column indexes
    // this mapping is given by the array of required indexes
    {
        IndexType numValues = haloJA.size();
        WriteAccess<IndexType> ja( haloJA );
        ReadAccess<IndexType> halo2global( halo.getRequiredIndexes() );

        for ( IndexType i = 0; i < numValues; i++ )
        {
            ja[i] = halo2global[ja[i]];
        }
    }
    LAMAArray<IndexType> outIA;
    LAMAArray<IndexType> outJA;
    LAMAArray<ValueType> outValues;
    IndexType numKeepDiagonals = 0;

    if ( attemptDiagonalProperty && localData.hasDiagonalProperty() )
    {
        numKeepDiagonals = std::min( localData.getNumRows(), localData.getNumColumns() );
        SCAI_LOG_INFO( logger, localData << ": has diagonal property, numKeepDiagonals = " << numKeepDiagonals );
    }

    // use static method of MatrixStorage
    StorageMethods<ValueType>::joinCSR( outIA, outJA, outValues, localIA, localJA, localValues, haloIA, haloJA,
                                        haloValues, numKeepDiagonals );
    // here mIA is size array, NOT offsets
    const IndexType numRows = outIA.size() - 1;
    const IndexType numColumns = colDist.getGlobalSize();
    const IndexType numValues = outJA.size();
    setCSRData( numRows, numColumns, numValues, outIA, outJA, outValues );
    check( "joined matrix storage" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::localize( const _MatrixStorage& globalData, const Distribution& rowDist )
{
    SCAI_REGION( "Storage.localize" )

    if ( rowDist.isReplicated() )
    {
        assign( globalData );
        return;
    }

    SCAI_ASSERT_EQUAL_ERROR( globalData.getNumRows(), rowDist.getGlobalSize() )
    const IndexType numColumns = globalData.getNumColumns();
    const IndexType localNumRows = rowDist.getLocalSize();
    LAMAArray<IndexType> globalIA;
    LAMAArray<IndexType> globalJA;
    LAMAArray<ValueType> globalValues;
    globalData.buildCSRData( globalIA, globalJA, globalValues );
    LAMAArray<IndexType> localIA;
    LAMAArray<IndexType> localJA;
    LAMAArray<ValueType> localValues;
    StorageMethods<ValueType>::localizeCSR( localIA, localJA, localValues, globalIA, globalJA, globalValues, rowDist );
    const IndexType localNumValues = localJA.size();
    setCSRData( localNumRows, numColumns, localNumValues, localIA, localJA, localValues );
    check( "localized matrix" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::replicate( const _MatrixStorage& localData, const Distribution& rowDist )
{
    SCAI_REGION( "Storage.replicate" )

    if ( rowDist.isReplicated() )
    {
        SCAI_LOG_INFO( logger, "replicate: assign due to replicated distribution" )
        assign( localData );
        return;
    }

    SCAI_ASSERT_EQUAL_ERROR( localData.getNumRows(), rowDist.getLocalSize() )
    const IndexType numColumns = localData.getNumColumns();
    const IndexType globalNumRows = rowDist.getGlobalSize();
    LAMAArray<IndexType> localIA;
    LAMAArray<IndexType> localJA;
    LAMAArray<ValueType> localValues;
    buildCSRData( localIA, localJA, localValues );
    LAMAArray<IndexType> globalIA;
    LAMAArray<IndexType> globalJA;
    LAMAArray<ValueType> globalValues;
    StorageMethods<ValueType>::replicateCSR( globalIA, globalJA, globalValues, localIA, localJA, localValues, rowDist );
    const IndexType globalNumValues = globalJA.size();
    setCSRData( globalNumRows, numColumns, globalNumValues, globalIA, globalJA, globalValues );
    check( "replicated matrix storage" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::splitHalo(
    MatrixStorage<ValueType>& localData,
    MatrixStorage<ValueType>& haloData,
    Halo& halo,
    const Distribution& colDist,
    const Distribution* rowDist ) const
{
    SCAI_REGION( "Storage.splitHalo" )
    SCAI_LOG_INFO( logger, *this << ": split according to column distribution " << colDist )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, colDist.getGlobalSize() )

    if ( getFormat() != Format::CSR )
    {
        LAMA_UNSUPPORTED( "splitHalo is not supported for " << getFormat() << ", converting to CSR!" );
    }

    if ( colDist.isReplicated() )
    {
        // if there is no column distribution, halo is not needed
        if ( rowDist )
        {
            localData.localize( *this, *rowDist );
        }
        else
        {
            localData.assign( *this );
        }

        haloData.allocate( mNumRows, 0 );
        halo = Halo(); // empty halo schedule
        return;
    }

    IndexType numRows = mNumRows;

    // check optional row distribution if specified

    if ( rowDist )
    {
        SCAI_LOG_INFO( logger, *this << ": split also localizes for " << *rowDist )
        SCAI_ASSERT_EQUAL_ERROR( mNumRows, rowDist->getGlobalSize() )
        numRows = rowDist->getLocalSize();
    }

    LAMAArray<IndexType> csrIA;
    LAMAArray<IndexType> csrJA;
    LAMAArray<ValueType> csrValues;
    buildCSRData( csrIA, csrJA, csrValues );
    SCAI_LOG_INFO( logger, *this << ": CSR data generated, #non-zeros = " << csrJA.size() )
    LAMAArray<IndexType> localIA;
    LAMAArray<IndexType> localJA;
    LAMAArray<ValueType> localValues;
    LAMAArray<IndexType> haloIA;
    LAMAArray<IndexType> haloJA;
    LAMAArray<ValueType> haloValues;
    StorageMethods<ValueType>::splitCSR( localIA, localJA, localValues, haloIA, haloJA, haloValues, csrIA, csrJA,
                                         csrValues, colDist, rowDist );
    SCAI_ASSERT_EQUAL_DEBUG( localIA.size(), numRows + 1 )
    SCAI_ASSERT_EQUAL_DEBUG( haloIA.size(), numRows + 1 )
    const IndexType haloNumValues = haloJA.size();
    const IndexType localNumValues = localJA.size();
    SCAI_LOG_INFO( logger,
                   *this << ": split into " << localNumValues << " local non-zeros " " and " << haloNumValues << " halo non-zeros" )
    const IndexType localNumColumns = colDist.getLocalSize();
    IndexType haloNumColumns; // will be available after remap
    // build the halo by the non-local indexes
    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );
    SCAI_LOG_INFO( logger, "build halo: " << halo )
    localData.setCSRData( numRows, localNumColumns, localNumValues, localIA, localJA, localValues );
    localData.check( "local part after split" );
    // halo data is expected to have many empty rows, so enable compressing with row indexes
    haloData.setCompressThreshold( 0.5 );
    haloData.setCSRData( numRows, haloNumColumns, haloNumValues, haloIA, haloJA, haloValues );
    haloData.check( "halo part after split" );
    SCAI_LOG_INFO( logger,
                   "Result of split: local storage = " << localData << ", halo storage = " << haloData << ", halo = " << halo )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::buildHalo( Halo& halo, const Distribution& colDist )
{
    SCAI_LOG_INFO( logger, *this << ": build halo according to column distribution " << colDist )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, colDist.getGlobalSize() )
    LAMAArray<IndexType> haloIA;
    LAMAArray<IndexType> haloJA; // global columns, all non-local
    LAMAArray<ValueType> haloValues;
    buildCSRData( haloIA, haloJA, haloValues );
    const IndexType haloNumValues = haloJA.size();
    IndexType haloNumColumns; // will be available after remap
    // build the halo by the non-local indexes
    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );
    setCSRData( mNumRows, haloNumColumns, haloNumValues, haloIA, haloJA, haloValues );
    check( "halo part after split" );
    SCAI_LOG_INFO( logger, "Result of buildHalo: " << "halo storage = " << *this << ", halo = " << halo )
}

/* --------------------------------------------------------------------------- */

void _MatrixStorage::scale( const ContextArray& )
{
    COMMON_THROWEXCEPTION( "scale of rows not supported yet, matrix = " << *this )
}

void _MatrixStorage::setDiagonal( const ContextArray& )
{
    COMMON_THROWEXCEPTION( "set Diagonal not suppported yet, matrix = " << *this )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::invert( const MatrixStorage<ValueType>& other )
{
    // there is only a feedback solution via dense matrices
    DenseStorage<ValueType> otherDense( other );
    otherDense.invert( otherDense );
    this->assign( otherDense );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::matrixTimesVector(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    LAMA_UNSUPPORTED( *this << ": no matrixTimesVector for this format available, take CSR" )
    CSRStorage<ValueType> tmp( *this );
    tmp.matrixTimesVector( result, alpha, x, beta, y );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::vectorTimesMatrix(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    LAMA_UNSUPPORTED( *this << ": no vectorTimesMatrix for this format available, take CSR" )
    CSRStorage<ValueType> tmp( *this );
    tmp.vectorTimesMatrix( result, alpha, x, beta, y );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::matrixTimesVectorN(
    LAMAArray<ValueType>& result,
    const IndexType n,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    LAMA_UNSUPPORTED(
        *this << ": no matrixTimesVectorN" << " ( denseStorage = anyStorage * denseStorage )" << " for this format available, take CSR" );
    CSRStorage<ValueType> tmp( *this );
    tmp.matrixTimesVectorN( result, n, alpha, x, beta, y );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* MatrixStorage<ValueType>::matrixTimesVectorAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger, *this << ": asynchronous matrixTimesVector by new thread" )
    // general default: asynchronous execution is done by a new thread
    void ( MatrixStorage::*pf )(
        LAMAArray<ValueType>&,
        const ValueType,
        const LAMAArray<ValueType>&,
        const ValueType,
        const LAMAArray<ValueType>& ) const
    = &MatrixStorage<ValueType>::matrixTimesVector;
    using scai::common::bind;
    using scai::common::ref;
    using scai::common::cref;
    return new TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* MatrixStorage<ValueType>::vectorTimesMatrixAsync(
    LAMAArray<ValueType>& result,
    const ValueType alpha,
    const LAMAArray<ValueType>& x,
    const ValueType beta,
    const LAMAArray<ValueType>& y ) const
{
    SCAI_LOG_INFO( logger, *this << ": asynchronous vectorTimesMatrix by new thread" )
    // general default: asynchronous execution is done by a new thread
    void ( MatrixStorage::*pf )(
        LAMAArray<ValueType>&,
        const ValueType,
        const LAMAArray<ValueType>&,
        const ValueType,
        const LAMAArray<ValueType>& ) const
    = &MatrixStorage<ValueType>::vectorTimesMatrix;
    using scai::common::bind;
    using scai::common::ref;
    using scai::common::cref;
    return new TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ) ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::jacobiIterate(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    LAMA_UNSUPPORTED( *this << ": no jacobiIterate for this format available, take CSR" )
    CSRStorage<ValueType> tmp( *this );
    tmp.jacobiIterate( solution, oldSolution, rhs, omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* MatrixStorage<ValueType>::jacobiIterateAsync(
    LAMAArray<ValueType>& solution,
    const LAMAArray<ValueType>& oldSolution,
    const LAMAArray<ValueType>& rhs,
    const ValueType omega ) const
{
    // general default: asynchronous execution is done by a new thread
    void ( MatrixStorage::*pf )(
        LAMAArray<ValueType>&,
        const LAMAArray<ValueType>&,
        const LAMAArray<ValueType>&,
        const ValueType ) const
    = &MatrixStorage<ValueType>::jacobiIterate;
    using scai::common::bind;
    using scai::common::cref;
    using scai::common::ref;
    return new TaskSyncToken( bind( pf, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_UNSUPPORTED( *this << ": jacobiIterateHalo for this format NOT available, take CSR" )
    CSRStorage<ValueType> tmpHalo( *this );
    // very inefficient as we just need the diagonal
    CSRStorage<ValueType> tmpLocal( localStorage );
    tmpHalo.jacobiIterateHalo( localSolution, tmpLocal, oldHaloSolution, omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::jacobiIterateHalo(
    LAMAArray<ValueType>& localSolution,
    const LAMAArray<ValueType>& localDiagonal,
    const LAMAArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    LAMA_UNSUPPORTED( *this << ": jacobiIterateHalo for this format NOT available, take CSR" )
    CSRStorage<ValueType> tmpHalo( *this );
    tmpHalo.jacobiIterateHalo( localSolution, localDiagonal, oldHaloSolution, omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::matrixTimesScalar( const ValueType alpha, const MatrixStorage<ValueType>& a )
{
    SCAI_LOG_INFO( logger, *this << " = alpha( " << alpha << " ) x " << a )

    if ( &a != this )
    {
        assign( a );
    }

    scale( alpha );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::matrixPlusMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const ValueType beta,
    const MatrixStorage<ValueType>& b )
{
    LAMA_UNSUPPORTED( *this << ": no matrixPlusMatrix for this format available, take CSR" )
    // TODO How can we make sure that CSRStorage really has overridden it, otherwise endless recursion here
    CSRStorage<ValueType> tmp( *this );
    tmp.check( "Temporary CSR storage for matrix addition" );
    tmp.matrixPlusMatrix( alpha, a, beta, b );
    assign( tmp );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::matrixTimesMatrix(
    const ValueType alpha,
    const MatrixStorage<ValueType>& a,
    const MatrixStorage<ValueType>& b,
    const ValueType beta,
    const MatrixStorage<ValueType>& y )
{
    LAMA_UNSUPPORTED( *this << ": no matrixTimesMatrix for this format available, take CSR" )
    // TODO How can we make sure that CSR really has overridden it, otherwise endless recursion here
    CSRStorage<ValueType> tmp( *this );
    tmp.check( "Temporary CSR storage for matrix multiplication" );
    tmp.matrixTimesMatrix( alpha, a, b, beta, y );
    assign( tmp );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
ValueType MatrixStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_ASSERT_EQUAL_ERROR( mNumRows, other.getNumRows() )
    SCAI_ASSERT_EQUAL_ERROR( mNumColumns, other.getNumColumns() )
    LAMA_UNSUPPORTED( *this << ": no maxDiffNorm for format " << getFormat() << " available, take Dense" )
    DenseStorage<ValueType> tmp( *this );
    return tmp.maxDiffNorm( other );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::exchangeHalo(
    const Halo& halo,
    const MatrixStorage<ValueType>& matrix,
    const Communicator& comm )
{
    IndexType numColumns = matrix.getNumColumns(); // remains unchanged
    // get the matrix data in CSR format
    LAMAArray<IndexType> sourceIA;
    LAMAArray<IndexType> sourceJA;
    LAMAArray<ValueType> sourceValues;
    matrix.buildCSRData( sourceIA, sourceJA, sourceValues );
    LAMAArray<IndexType> targetIA;
    LAMAArray<IndexType> targetJA;
    LAMAArray<ValueType> targetValues;
    StorageMethods<ValueType>::exchangeHaloCSR( targetIA, targetJA, targetValues, sourceIA, sourceJA, sourceValues,
            halo, comm );
    const IndexType targetNumRows = targetIA.size() - 1;
    const IndexType targetNumValues = targetJA.size();
    setCSRData( targetNumRows, numColumns, targetNumValues, targetIA, targetJA, targetValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::redistribute( const _MatrixStorage& other, const Redistributor& redistributor )

{
    if ( other.getFormat() == Format::CSR && other.getValueType() == getValueType() )
    {
        // This special case avoids unnecssary CSR conversions
        const CSRStorage<ValueType>* otherCSR = dynamic_cast<const CSRStorage<ValueType>*>( &other );
        SCAI_ASSERT_DEBUG( otherCSR, "serious cast error" );
        redistributeCSR( *otherCSR, redistributor );
        return;
    }

    SCAI_REGION( "Storage.redistribute" )
    // For the redistribution we use the CSR format on both sides
    const Distribution& sourceDistribution = *redistributor.getSourceDistributionPtr();
    const Distribution& targetDistribution = *redistributor.getTargetDistributionPtr();
    SCAI_LOG_INFO( logger, other << ": redistribute rows via " << redistributor )
    bool sameDist = false;

    // check for same distribution, either equal or both replicated

    if ( sourceDistribution.isReplicated() && targetDistribution.isReplicated() )
    {
        sameDist = true;
    }
    else if ( &sourceDistribution == &targetDistribution )
    {
        sameDist = true;
    }

    if ( sameDist )
    {
        SCAI_LOG_INFO( logger, "redistributor with same source/target distribution" )
        assign( other );
        return; // so we are done
    }

    const IndexType numColumns = other.getNumColumns(); // does not change

    // check that source distribution fits with storage
    SCAI_ASSERT_EQUAL_ERROR( other.getNumRows(), sourceDistribution.getLocalSize() )
    // get the matrix data from other in CSR format
    LAMAArray<IndexType> sourceIA;

    LAMAArray<IndexType> sourceJA;

    LAMAArray<ValueType> sourceValues;

    other.buildCSRData( sourceIA, sourceJA, sourceValues );

    LAMAArray<IndexType> targetIA;

    LAMAArray<IndexType> targetJA;

    LAMAArray<ValueType> targetValues;

    StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues, sourceIA, sourceJA, sourceValues,
            redistributor );

    const IndexType targetNumRows = targetIA.size() - 1;

    const IndexType targetNumValues = targetJA.size();

    setCSRData( targetNumRows, numColumns, targetNumValues, targetIA, targetJA, targetValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::redistributeCSR( const CSRStorage<ValueType>& other, const Redistributor& redistributor )
{
    SCAI_REGION( "Storage.redistributeCSR" )
    const Distribution& sourceDistribution = *redistributor.getSourceDistributionPtr();
    const Distribution& targetDistribution = *redistributor.getTargetDistributionPtr();
    SCAI_LOG_INFO( logger, other << ": redistribute rows via " << redistributor )
    bool sameDist = false;

    // check for same distribution, either equal or both replicated

    if ( sourceDistribution.isReplicated() && targetDistribution.isReplicated() )
    {
        sameDist = true;
    }
    else if ( &sourceDistribution == &targetDistribution )
    {
        sameDist = true;
    }

    if ( sameDist )
    {
        SCAI_LOG_INFO( logger, "redistributor with same source/target distribution" )
        assign( other );
        return; // so we are done
    }

    const IndexType numColumns = other.getNumColumns(); // does not change

    // check that source distribution fits with storage
    SCAI_ASSERT_EQUAL_ERROR( other.getNumRows(), sourceDistribution.getLocalSize() )
    // it is not necessary to convert the other storage to CSR
    LAMAArray<IndexType> targetIA;

    LAMAArray<IndexType> targetJA;

    LAMAArray<ValueType> targetValues;

    StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues, other.getIA(), other.getJA(),
            other.getValues(), redistributor );

    const IndexType targetNumRows = targetIA.size() - 1;

    const IndexType targetNumValues = targetJA.size();

    setCSRData( targetNumRows, numColumns, targetNumValues, targetIA, targetJA, targetValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
template<typename OtherValueType>
void MatrixStorage<ValueType>::setRawDenseData(
    const IndexType numRows,
    const IndexType numColumns,
    const OtherValueType values[],
    const ValueType epsilon )
{
    SCAI_ASSERT_ERROR( epsilon >= 0.0, "epsilon = " << epsilon << ", must not be negative" )
    mEpsilon = epsilon;
    // wrap all the data in a dense storage and make just an assign
    SCAI_LOG_INFO( logger, "set dense storage " << numRows << " x " << numColumns )
    LAMAArrayRef<OtherValueType> data( numRows * numColumns, values );
    SCAI_LOG_INFO( logger, "use LAMA array ref: " << data << ", size = " << data.size() )
    DenseStorageView<OtherValueType> denseStorage( data, numRows, numColumns );
    assign( denseStorage ); // will internally use the value epsilon
    SCAI_LOG_INFO( logger, *this << ": have set dense data " << numRows << " x " << numColumns )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::setDenseData(
    const IndexType numRows,
    const IndexType numColumns,
    const ContextArray& values,
    const ValueType epsilon )
{
    mEpsilon = epsilon;
    // const_cast required, is safe as we will create a const DenseStorageView
    ContextArray& mValues = const_cast<ContextArray&>( values );

    switch ( values.getValueType() )
    {

#define LAMA_DENSE_ASSIGN( z, I, _ )                                                             \
case common::scalar::SCALAR_ARITHMETIC_TYPE##I :                                                 \
{                                                                                                \
    LAMAArray<ARITHMETIC_TYPE##I>& typedValues =                                                 \
            dynamic_cast<LAMAArray<ARITHMETIC_TYPE##I>&>( mValues );                             \
    const DenseStorageView<ARITHMETIC_TYPE##I> denseStorage( typedValues, numRows, numColumns ); \
    ARITHMETIC_TYPE##I tmpEpsilon = static_cast<ARITHMETIC_TYPE##I>( epsilon );                  \
    denseStorage.swapEpsilon( tmpEpsilon );                                                      \
    assign( denseStorage );                                                                      \
    break;                                                                                       \
}

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_DENSE_ASSIGN, _ )

#undef LAMA_DENSE_ASSIGN

        default        :
            COMMON_THROWEXCEPTION( "Unsupported type for setting dense data: " << values.getValueType() )
    }
}

/* ========================================================================= */
/*       File I/O                                                            */
/* ========================================================================= */

template<typename ValueType>
void MatrixStorage<ValueType>::writeToFile(
    const std::string& fileName,
    const File::FileType fileType,
    const File::DataType dataType,
    const File::IndexDataType indexDataTypeIA,
    const File::IndexDataType indexDataTypeJA ) const
{
    writeToFile( 1, 0, fileName, fileType, dataType, indexDataTypeIA, indexDataTypeJA );
}

template<typename ValueType>
void MatrixStorage<ValueType>::writeToFile(
    const PartitionId size,
    const PartitionId rank,
    const std::string& fileName,
    const File::FileType fileType,
    const File::DataType dataType,
    const File::IndexDataType indexDataTypeIA,
    const File::IndexDataType indexDataTypeJA ) const
{
    LAMAArray<IndexType> csrIA;
    LAMAArray<IndexType> csrJA;
    LAMAArray<ValueType> csrValues;
// TODO Do not build CSR if this matrix is CSR storage
    buildCSRData( csrIA, csrJA, csrValues );
    StorageIO<ValueType>::writeCSRToFile( size, rank, csrIA, mNumColumns, csrJA, csrValues, fileName, fileType,
                                          dataType, indexDataTypeIA, indexDataTypeJA );
}

/*****************************************************************************/

template<typename ValueType>
void MatrixStorage<ValueType>::readFromFile( const std::string& fileName )
{
    SCAI_LOG_INFO( logger, "MatrixStorage<" << getValueType() << ">::readFromFile( " << fileName << ")" )
    SCAI_REGION( "Storage.readFromFile" )
    IndexType numColumns;
    IndexType numRows;
    IndexType numValues;
    LAMAArray<IndexType> csrIA;
    LAMAArray<IndexType> csrJA;
    LAMAArray<ValueType> csrValues;
    StorageIO<ValueType>::readCSRFromFile( csrIA, numColumns, csrJA, csrValues, fileName );
    numRows = csrIA.size() - 1;
    numValues = csrJA.size();
    SCAI_LOG_INFO( logger,
                   "read CSR storage <" << getValueType() << "> : " << numRows << " x " << numColumns << ", #values = " << numValues )
    setCSRData( numRows, numColumns, numValues, csrIA, csrJA, csrValues );
    check( "read matrix" );
}

/*****************************************************************************/

void _MatrixStorage::buildCSRGraph(
    IndexType* adjIA,
    IndexType* adjJA,
    IndexType* vwgt,
    CommunicatorPtr comm,
    const IndexType* globalRowIndexes /* = NULL */,
    IndexType* vtxdist /* = NULL */ ) const
{
    IndexType numLocalRows = mNumRows;

    if ( vtxdist != NULL ) // parallel graph
    {
        const PartitionId MASTER = 0;
        IndexType parts = comm->getSize();
// Is this valid ?
// SCAI_ASSERT_ERROR( getDistribution().getNumPartitions() == parts,
//              "mismatch number of partitions and communicator size" );
        std::vector<IndexType> localNumRows( parts );
        comm->gather( vtxdist, 1, MASTER, &numLocalRows );
        comm->bcast( vtxdist, parts, MASTER );
        vtxdist[parts] = OpenMPCSRUtils::scan( vtxdist, parts );
    }

    LAMAArray<IndexType> csrIA;
    LAMAArray<IndexType> csrJA;
    LAMAArray<float> csrValues;
    buildCSRData( csrIA, csrJA, csrValues );
    LAMAArray<IndexType> rowSizes;
    WriteOnlyAccess<IndexType> sizes( rowSizes, mNumRows );
    ReadAccess<IndexType> ia( csrIA );
    OpenMPCSRUtils::offsets2sizes( sizes.get(), ia.get(), mNumRows );
    ReadAccess<IndexType> ja( csrJA );
    IndexType offset = 0; // runs through JA
    IndexType newOffset = 0; // runs through adjJA

    for ( IndexType i = 0; i < numLocalRows; i++ )
    {
        IndexType rowIndex = i;

        if ( globalRowIndexes != NULL )
        {
            rowIndex = globalRowIndexes[i];
        }

        adjIA[i] = newOffset;
        vwgt[i] = sizes[i];

        for ( IndexType jj = 0; jj < sizes[i]; jj++ )
        {
            if ( rowIndex != ja[offset] ) // skip diagonal element
            {
                adjJA[newOffset++] = ja[offset];
            }

            offset++;
        }
    }

    adjIA[numLocalRows] = newOffset;
}

template<typename ValueType>
bool MatrixStorage<ValueType>::checkSymmetry() const
{
// check symmetry of matrix
    IndexType n = mNumRows;

    if ( n != mNumColumns )
    {
        return false;
    }

    for ( IndexType i = 0; i < n; ++i )
    {
        for ( IndexType j = 0; j < i; ++j )
        {
            if ( getValue( i, j ) != getValue( j, i ) )
            {
                return false;
            }
        }
    }

    return true;
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_MATRIX_STORAGE2_INSTANTIATE(z, J, TYPE )      \
    template COMMON_DLL_IMPORTEXPORT                             \
    void MatrixStorage<TYPE>::setRawDenseData(                 \
            const IndexType numRows,                               \
            const IndexType numColumns,                            \
            const ARITHMETIC_TYPE##J values[],                     \
            const TYPE );

#define LAMA_MATRIX_STORAGE_INSTANTIATE(z, I, _)                          \
    \
    template class COMMON_DLL_IMPORTEXPORT MatrixStorage<ARITHMETIC_TYPE##I> ;  \
    \
    BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT,                                 \
                     LAMA_MATRIX_STORAGE2_INSTANTIATE,                    \
                     ARITHMETIC_TYPE##I )                                 \
     

BOOST_PP_REPEAT( ARITHMETIC_TYPE_CNT, LAMA_MATRIX_STORAGE_INSTANTIATE, _ )

#undef LAMA_MATRIX_STORAGE_INSTANTIATE
#undef LAMA_MATRIX_STORAGE2_INSTANTIATE

} /* end namespace lama */

} /* end namespace scai */

std::ostream& operator<<( std::ostream& stream, const scai::lama::MatrixStorageFormat storageFormat )
{
    stream << scai::lama::format2Str( storageFormat );
    return stream;
}
