/**
 * @file MatrixStorage.cpp
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
 * @brief Implementation of methods for common base class of all matrix storage formats.
 * @author Thomas Brandes
 * @date 27.04.2011
 */

// hpp
#include <scai/lama/storage/MatrixStorage.hpp>

// local library
#include <scai/lama/storage/CSRStorage.hpp>
#include <scai/lama/storage/DenseStorage.hpp>
#include <scai/lama/storage/StorageMethods.hpp>

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/Redistributor.hpp>
#include <scai/dmemo/Halo.hpp>

#include <scai/lama/io/FileIO.hpp>

// internal scai libraries
#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/UtilKernelTrait.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/openmp/OpenMPUtils.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/tracing.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/unsupported.hpp>
#include <scai/common/macros/instantiate.hpp>
#include <scai/common/macros/loop.hpp>

#include <functional>

namespace scai
{

using namespace hmemo;
using namespace dmemo;

using tasking::SyncToken;
using tasking::TaskSyncToken;

using utilskernel::LAMAKernel;
using utilskernel::UtilKernelTrait;
using utilskernel::OpenMPUtils;

using sparsekernel::CSRKernelTrait;
using sparsekernel::OpenMPCSRUtils;

using common::BinaryOp;

namespace lama
{

/* --------------------------------------------------------------------------- */

template<typename ValueType>
MatrixStorage<ValueType>::MatrixStorage( const IndexType numRows, const IndexType numColumns, ContextPtr ctx ) : 

    _MatrixStorage( numRows, numColumns, ctx )

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
MatrixStorage<ValueType>::MatrixStorage( const MatrixStorage<ValueType>& other ) : 

    _MatrixStorage( other )
{
}

template<typename ValueType>
MatrixStorage<ValueType>::MatrixStorage( MatrixStorage<ValueType>&& other ) noexcept :

    _MatrixStorage( std::move( other ) )

{
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
common::ScalarType MatrixStorage<ValueType>::getValueType() const
{
    return common::getScalarType<ValueType>();
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::convertCSR2CSC(
    HArray<IndexType>& colIA,
    HArray<IndexType>& colJA,
    HArray<ValueType>& colValues,
    const IndexType numColumns,
    const HArray<IndexType>& rowIA,
    const HArray<IndexType>& rowJA,
    const HArray<ValueType>& rowValues,
    const ContextPtr preferredLoc )
{
    const IndexType numRows = rowIA.size() - 1;
    const IndexType numValues = rowJA.size();
    SCAI_ASSERT_EQUAL_DEBUG( rowJA.size(), rowValues.size() )
    static LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = preferredLoc;
    convertCSR2CSC.getSupportedContext( loc );
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
    convertCSR2CSC[loc]( cIA.get(), cJA.get(), cValues.get(),  // output args
                         rIA.get(), rJA.get(), rValues.get(), numRows, numColumns, numValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::buildCSCData(
    HArray<IndexType>& colIA,
    HArray<IndexType>& colJA,
    HArray<ValueType>& colValues ) const
{
    HArray<IndexType> rowIA;
    HArray<IndexType> rowJA;
    HArray<ValueType> rowValues;
    buildCSRData( rowIA, rowJA, rowValues );
    ContextPtr loc = Context::getHostPtr();
    convertCSR2CSC( colIA, colJA, colValues, getNumColumns(), rowIA, rowJA, rowValues, loc );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::getFirstColumnIndexes( hmemo::HArray<IndexType>& colIndexes ) const
{
    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    // ToDo: csrIA[i] == csrIA[i+1], no entry in row at all
    // ToDo: csrIA[numRows-1] == numValues possible, out of range addressing

    // gather: colIndexes[i] = csrJA[ csrIA[i] ]

    if ( getNumRows() > 0 )
    {
        SCAI_ASSERT_LT_ERROR( csrIA[getNumRows() - 1], csrJA.size(), "last row without any entry" )
    }

    static LAMAKernel<UtilKernelTrait::setGather<IndexType, IndexType> > setGather;

    ContextPtr loc = getContextPtr();
    setGather.getSupportedContext( loc );

    WriteOnlyAccess<IndexType> wColIndexes( colIndexes, loc, getNumRows() );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> ja( csrJA, loc );
    ReadAccess<IndexType> ia( csrIA, loc );
    setGather[loc] ( wColIndexes.get(), ja.get(), ia.get(), BinaryOp::COPY, getNumRows() );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::assignDummy( const _MatrixStorage& other )
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

    if ( common::typeSize( other.getValueType() ) < common::typeSize( getValueType() ) )
    {
        other.copyTo( *this );
        return;
    }

    SCAI_LOG_INFO( logger, *this << ": (default) assign ( " << other << " )" )

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    other.buildCSRData( csrIA, csrJA, csrValues );

    IndexType numRows = other.getNumRows();
    IndexType numColumns = other.getNumColumns();

    SCAI_LOG_DEBUG( logger, "build CSR data " << numRows << " x " << numColumns << ", #nnz = " << csrJA.size() )

    setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );

    SCAI_LOG_DEBUG( logger, "now assigned: " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::copyTo( _MatrixStorage& other ) const
{
    SCAI_REGION( "Storage.copyTo" )

    // convert my matrix data to CSR and set it for the other matrix 

    SCAI_LOG_INFO( logger, *this << ": (default) copyTo ( " << other << " )" )

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    SCAI_ASSERT_EQ_DEBUG( csrIA.size(), getNumRows() + 1, "csr: offset array ia has illegal size" )
    SCAI_ASSERT_EQ_DEBUG( csrValues.size(), csrJA.size(), 
                          "csr: ja and values must have same number of entries for non-zero values" )

    SCAI_LOG_DEBUG( logger, "build CSR data " << getNumRows() << " x " << getNumColumns() << ", #nnz = " << csrJA.size() )

    other.setCSRData( getNumRows(), getNumColumns(), csrIA, csrJA, csrValues );

    SCAI_LOG_INFO( logger, "now assigned: " << *this )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::copyBlockTo( _MatrixStorage& other, const IndexType first, const IndexType n ) const
{
    using namespace utilskernel;

    SCAI_ASSERT_LE( first, first + n, "illegal range" )
    SCAI_ASSERT_VALID_INDEX( first, getNumRows(), "first row out of range" )
    SCAI_ASSERT_VALID_INDEX( first + n - 1, getNumRows(), "last row out of range" );

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    ContextPtr loc = this->getContextPtr();

    SCAI_LOG_INFO( logger, "copyBlockTo : first = " << first << ", n = " << n << ", from this : " << *this )

    // copy out the corresponding sections, ia needs a shifting to zero

    HArray<IndexType> blockIA( n + 1 );
    HArrayUtils::setArraySection( blockIA, 0, 1, csrIA, first, 1, n +  1, BinaryOp::COPY, loc );

    IndexType offset = blockIA[0];  // gives shifting, as blockIA[0] must be 0
    HArrayUtils::compute( blockIA, blockIA, BinaryOp::SUB, offset, loc );

    IndexType numBlockValues = blockIA[n];

    SCAI_LOG_DEBUG( logger, "offset = " << offset << ", #nnz = " << numBlockValues );

    HArray<IndexType> blockJA( numBlockValues );
    HArray<ValueType> blockValues( numBlockValues );

    HArrayUtils::setArraySection( blockJA, 0, 1, csrJA, offset, 1, numBlockValues, BinaryOp::COPY, loc );
    HArrayUtils::setArraySection( blockValues, 0, 1, csrValues, offset, 1, numBlockValues, BinaryOp::COPY, loc );

    other.setCSRData( n, getNumColumns(), blockIA, blockJA, blockValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::assignTranspose( const MatrixStorage<ValueType>& other )
{
    SCAI_REGION( "Storage.assignTranspose" )
    SCAI_LOG_INFO( logger, *this << ": assignTranspose " << other )
    HArray<IndexType> cscIA;
    HArray<IndexType> cscJA;
    HArray<ValueType> cscValues;
    other.buildCSCData( cscIA, cscJA, cscValues );
    // Compressed sparse column data can be used directly to generate the transposed matrix
    // by interpretation as CSR data.
    setCSRData( other.getNumColumns(), other.getNumRows(), cscIA, cscJA, cscValues );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::joinRows(
    HArray<IndexType>& outSizes,
    HArray<IndexType>& outJA,
    HArray<ValueType>& outValues,
    const IndexType numLocalRows,
    const HArray<IndexType>& rowIndexes,
    const HArray<IndexType>& inSizes,
    const HArray<IndexType>& inJA,
    const HArray<ValueType>& inValues )
{
    SCAI_REGION( "Storage.joinRows" )
    SCAI_LOG_INFO( logger, "join " << numLocalRows << " rows " )
    {
        WriteOnlyAccess<IndexType> sizes( outSizes, numLocalRows );
        ReadAccess<IndexType> rowSizes( inSizes );
        ReadAccess<IndexType> indexes( rowIndexes );

        // initialize counters (Attention: sizes.size() != rowSizes.size())

        for ( IndexType i = 0; i < sizes.size(); ++i )
        {
            sizes[i] = 0;
        }

        // count elements for each row
        for ( IndexType i = 0; i < rowSizes.size(); ++i )
        {
            sizes[indexes[i]] += rowSizes[i];
        }
    }
    // generate offset array for insertion
    HArray<IndexType> IA;
    {
        WriteOnlyAccess<IndexType> offsets( IA, numLocalRows + 1 );
        ReadAccess<IndexType> sizes( outSizes );
        OpenMPUtils::set( offsets.get(), sizes.get(), numLocalRows, BinaryOp::COPY );
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

    for ( IndexType i = 0; i < rowSizes.size(); ++i )
    {
        IndexType currentRow = indexes[i];

        for ( IndexType ii = 0; ii < rowSizes[i]; ++ii )
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
MatrixStorage<ValueType>& MatrixStorage<ValueType>::operator=( const MatrixStorage<ValueType>& other )
{
    assign( other ); // assign can deal with all kind of storage formats/types
    return *this;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::moveImpl( MatrixStorage<ValueType>&& other )
{
    _MatrixStorage::moveImpl( std::move( other ) );
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
    HArray<IndexType> localIA;
    HArray<IndexType> localJA;
    HArray<ValueType> localValues;
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
    HArray<IndexType> haloIA;
    HArray<IndexType> haloJA;
    HArray<ValueType> haloValues;
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
    HArray<IndexType> outIA;
    HArray<IndexType> outJA;
    HArray<ValueType> outValues;
    IndexType numKeepDiagonals = 0;

    if ( attemptDiagonalProperty && localData.hasDiagonalProperty() )
    {
        numKeepDiagonals = common::Math::min( localData.getNumRows(), localData.getNumColumns() );
        SCAI_LOG_INFO( logger, localData << ": has diagonal property, numKeepDiagonals = " << numKeepDiagonals );
    }

    // use static method of MatrixStorage
    StorageMethods<ValueType>::joinCSR( outIA, outJA, outValues, localIA, localJA, localValues, haloIA, haloJA,
                                        haloValues, numKeepDiagonals );
    // here mIA is size array, NOT offsets
    const IndexType numRows = outIA.size() - 1;
    const IndexType numColumns = colDist.getGlobalSize();
    setCSRData( numRows, numColumns, outIA, outJA, outValues );
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
    HArray<IndexType> globalIA;
    HArray<IndexType> globalJA;
    HArray<ValueType> globalValues;
    globalData.buildCSRData( globalIA, globalJA, globalValues );
    HArray<IndexType> localIA;
    HArray<IndexType> localJA;
    HArray<ValueType> localValues;
    StorageMethods<ValueType>::localizeCSR( localIA, localJA, localValues, globalIA, globalJA, globalValues, rowDist );
    setCSRData( localNumRows, numColumns, localIA, localJA, localValues );
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
    HArray<IndexType> localIA;
    HArray<IndexType> localJA;
    HArray<ValueType> localValues;
    buildCSRData( localIA, localJA, localValues );
    HArray<IndexType> globalIA;
    HArray<IndexType> globalJA;
    HArray<ValueType> globalValues;
    StorageMethods<ValueType>::replicateCSR( globalIA, globalJA, globalValues, localIA, localJA, localValues, rowDist );
    setCSRData( globalNumRows, numColumns, globalIA, globalJA, globalValues );
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
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), colDist.getGlobalSize() )

    if ( getFormat() != Format::CSR )
    {
        SCAI_UNSUPPORTED( "splitHalo is not supported for " << getFormat() << ", converting to CSR!" );
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

        haloData.allocate( getNumRows(), 0 );
        halo = Halo(); // empty halo schedule
        return;
    }

    IndexType numRows = getNumRows();

    // check optional row distribution if specified

    if ( rowDist )
    {
        SCAI_LOG_INFO( logger, *this << ": split also localizes for " << *rowDist )
        SCAI_ASSERT_EQUAL_ERROR( getNumRows(), rowDist->getGlobalSize() )
        numRows = rowDist->getLocalSize();
    }

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;
    buildCSRData( csrIA, csrJA, csrValues );
    SCAI_LOG_INFO( logger, *this << ": CSR data generated, #non-zeros = " << csrJA.size() )
    HArray<IndexType> localIA;
    HArray<IndexType> localJA;
    HArray<ValueType> localValues;
    HArray<IndexType> haloIA;
    HArray<IndexType> haloJA;
    HArray<ValueType> haloValues;
    StorageMethods<ValueType>::splitCSR( localIA, localJA, localValues, haloIA, haloJA, haloValues, csrIA, csrJA,
                                         csrValues, colDist, rowDist );
    SCAI_ASSERT_EQUAL_DEBUG( localIA.size(), numRows + 1 )
    SCAI_ASSERT_EQUAL_DEBUG( haloIA.size(), numRows + 1 )
    SCAI_LOG_INFO( logger,
                   *this << ": split into " << localJA.size() << " local non-zeros " " and " << haloJA.size() << " halo non-zeros" )
    const IndexType localNumColumns = colDist.getLocalSize();
    IndexType haloNumColumns; // will be available after remap
    // build the halo by the non-local indexes
    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );
    SCAI_LOG_INFO( logger, "build halo: " << halo )
    localData.setCSRData( numRows, localNumColumns, localIA, localJA, localValues );
    localData.check( "local part after split" );
    // halo data is expected to have many empty rows, so enable compressing with row indexes
    haloData.setCompressThreshold( 0.5 );
    haloData.setCSRData( numRows, haloNumColumns, haloIA, haloJA, haloValues );
    haloData.check( "halo part after split" );
    SCAI_LOG_INFO( logger,
                   "Result of split: local storage = " << localData << ", halo storage = " << haloData << ", halo = " << halo )
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::buildHalo( Halo& halo, const Distribution& colDist )
{
    SCAI_LOG_INFO( logger, *this << ": build halo according to column distribution " << colDist )
    SCAI_ASSERT_EQUAL_ERROR( getNumColumns(), colDist.getGlobalSize() )
    HArray<IndexType> haloIA;
    HArray<IndexType> haloJA; // global columns, all non-local
    HArray<ValueType> haloValues;
    buildCSRData( haloIA, haloJA, haloValues );
    IndexType haloNumColumns; // will be available after remap
    // build the halo by the non-local indexes
    _StorageMethods::buildHalo( halo, haloJA, haloNumColumns, colDist );
    setCSRData( getNumRows(), haloNumColumns, haloIA, haloJA, haloValues );
    check( "halo part after split" );
    SCAI_LOG_INFO( logger, "Result of buildHalo: " << "halo storage = " << *this << ", halo = " << halo )
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::compress( const RealType<ValueType> eps, bool keepDiagonal )
{
    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    CSRStorage<ValueType>::compress( csrIA, csrJA, csrValues, keepDiagonal, eps, getContextPtr() );

    setCSRData( getNumRows(), getNumColumns(), csrIA, csrJA, csrValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::globalizeHaloIndexes( const dmemo::Halo& halo, const IndexType globalNumColumns )
{
    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    halo.halo2Global( csrJA );

    setCSRData( getNumRows(), globalNumColumns, csrIA, csrJA, csrValues );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::invert( const MatrixStorage<ValueType>& other )
{
    SCAI_ASSERT_NE_ERROR( other.getFormat(), Format::DENSE, "default implementation relies on override for DenseStorage" )

    // there is only a fallback solution via dense matrices

    auto otherDense = convert<DenseStorage<ValueType>>( other );
    otherDense.invert( otherDense );
    this->assign( otherDense );
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::reduce(
    hmemo::HArray<ValueType>& array, 
    const IndexType dim, 
    const common::BinaryOp reduceOp,
    const common::UnaryOp elemOp )
{
    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    hmemo::ReadAccess<ValueType> values( csrValues );
    hmemo::ReadAccess<IndexType> ia( csrIA );
    hmemo::ReadAccess<IndexType> ja( csrJA );

    hmemo::WriteAccess<ValueType> a( array );

    OpenMPCSRUtils::reduce( a.get(), ia.get(), ja.get(), values.get(), getNumRows(), dim, reduceOp, elemOp );

    /*
    if ( dim == 0 )
    {
        SCAI_ASSERT_EQ_ERROR( array.size(), getNumRows(), "size mismatch" );

        for ( IndexType i = 0; i < getNumRows(); ++i )
        {
            // a[i] = 0;

            for ( IndexType jj = ia[i]; jj < ia[i+1]; ++jj )
            {
                ValueType v = values[jj];
                v = applyUnaryOp( elemOp, v );
                a[i] += v;
            }
        }
    }
    else if ( dim == 1 )
    {
        SCAI_ASSERT_EQ_ERROR( array.size(), getNumColumns(), "size mismatch" );

        // for ( IndexType j = 0; j < getNumColumns(); ++j )
        // {
            // a[j] = 0;
        // }

        for ( IndexType i = 0; i < getNumRows(); ++i )
        {
            for ( IndexType jj = ia[i]; jj < ia[i+1]; ++jj )
            {
                IndexType j = ja[jj];
                ValueType v = values[jj];
                v = applyUnaryOp( elemOp, v );
                a[j] += v;
            }
        }
    }
    */
}

/* ------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::matrixTimesVectorN(
    HArray<ValueType>& result,
    const IndexType n,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y ) const
{
    SCAI_UNSUPPORTED(
        *this << ": no matrixTimesVectorN" << " ( denseStorage = anyStorage * denseStorage )" << " for this format available, take CSR" );

    auto tmp = convert<CSRStorage<ValueType>>( *this );
    tmp.matrixTimesVectorN( result, n, alpha, x, beta, y );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* MatrixStorage<ValueType>::matrixTimesVectorAsync(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const common::MatrixOp op ) const
{
    SCAI_LOG_INFO( logger, *this << ": asynchronous matrixTimesVector by new thread" )
    // general default: asynchronous execution is done by a new thread
    void ( MatrixStorage::*pf )(
        HArray<ValueType>&,
        const ValueType,
        const HArray<ValueType>&,
        const ValueType,
        const HArray<ValueType>&,
        const common::MatrixOp ) const
    = &MatrixStorage<ValueType>::matrixTimesVector;
    using std::bind;
    using std::ref;
    using std::cref;
    return new TaskSyncToken( bind( pf, this, ref( result ), alpha, cref( x ), beta, cref( y ), op ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::jacobiIterate(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    SCAI_UNSUPPORTED( *this << ": no jacobiIterate for this format available, take CSR" )
    auto tmp = convert<CSRStorage<ValueType>>( * this );
    tmp.jacobiIterate( solution, oldSolution, rhs, omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* MatrixStorage<ValueType>::jacobiIterateAsync(
    HArray<ValueType>& solution,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const ValueType omega ) const
{
    // general default: asynchronous execution is done by a new thread
    void ( MatrixStorage::*pf )(
        HArray<ValueType>&,
        const HArray<ValueType>&,
        const HArray<ValueType>&,
        const ValueType ) const
    = &MatrixStorage<ValueType>::jacobiIterate;
    using std::bind;
    using std::cref;
    using std::ref;
    return new TaskSyncToken( bind( pf, this, ref( solution ), cref( oldSolution ), cref( rhs ), omega ) );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const MatrixStorage<ValueType>& localStorage,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_UNSUPPORTED( *this << ": jacobiIterateHalo for this format NOT available, take CSR" )

    // very inefficient as we just need the diagonal

    auto csrHalo  = convert<CSRStorage<ValueType>>( *this );
   
    // we need the diagonal from the local storage

    HArray<ValueType> diagValues;
    localStorage.getDiagonal( diagValues );

    auto csrLocal = diagonal<CSRStorage<ValueType>>( diagValues );

    csrHalo.jacobiIterateHalo( localSolution, csrLocal, oldHaloSolution, omega );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::jacobiIterateHalo(
    HArray<ValueType>& localSolution,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const ValueType omega ) const
{
    SCAI_UNSUPPORTED( *this << ": jacobiIterateHalo for this format NOT available, take CSR" )

    auto csrHalo  = convert<CSRStorage<ValueType>>( *this );
    csrHalo.jacobiIterateHalo( localSolution, localDiagonal, oldHaloSolution, omega );
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
    SCAI_UNSUPPORTED( *this << ": no matrixPlusMatrix for this format available, take CSR" )

    // Make sure that CSRStorage really has overridden it, otherwise endless recursion here

    SCAI_ASSERT_NE_ERROR( getFormat(), Format::CSR, "default implementation has not been overridden by CSR" )

    CSRStorage<ValueType> csr;
    csr.matrixPlusMatrix( alpha, a, beta, b );
    assign( csr );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void MatrixStorage<ValueType>::binaryOp(
    const MatrixStorage<ValueType>& a,
    const common::BinaryOp op,
    const MatrixStorage<ValueType>& b )
{
    // Make sure that CSRStorage really has overridden it, otherwise endless recursion here

    SCAI_ASSERT_NE_ERROR( getFormat(), Format::CSR, "default implementation has not been overridden by CSR" )

    CSRStorage<ValueType> csr;
    csr.binaryOp( a, op, b );
    assign( csr );
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
    SCAI_UNSUPPORTED( *this << ": no matrixTimesMatrix for this format available, take CSR" )

    // Make sure that CSRStorage really has overridden it, otherwise endless recursion here

    SCAI_ASSERT_NE_ERROR( getFormat(), Format::CSR, "default implementation has not been overridden by CSR" )

    auto csr  = convert<CSRStorage<ValueType>>( *this );
    csr.matrixTimesMatrix( alpha, a, b, beta, y );
    assign( csr );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> MatrixStorage<ValueType>::maxDiffNorm( const MatrixStorage<ValueType>& other ) const
{
    SCAI_ASSERT_EQ_ERROR( getNumRows(), other.getNumRows(), "row size mismatch for maxDiffNorm" )
    SCAI_ASSERT_EQ_ERROR( getNumColumns(), other.getNumColumns(), "col size mismatch for maxDiffNorm" )

    SCAI_UNSUPPORTED( *this << ": no maxDiffNorm for format " << getFormat() << " available, take Dense" )

    SCAI_ASSERT_NE_ERROR( getFormat(), Format::DENSE, "default implementation has not been overridden for DENSE" )

    auto dense = convert<DenseStorage<ValueType>>( *this );
    return dense.maxDiffNorm( other );
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
    HArray<IndexType> sourceIA;
    HArray<IndexType> sourceJA;
    HArray<ValueType> sourceValues;
    matrix.buildCSRData( sourceIA, sourceJA, sourceValues );
    HArray<IndexType> targetIA;
    HArray<IndexType> targetJA;
    HArray<ValueType> targetValues;
    StorageMethods<ValueType>::exchangeHaloCSR( targetIA, targetJA, targetValues, sourceIA, sourceJA, sourceValues,
            halo, comm );
    const IndexType targetNumRows = targetIA.size() - 1;
    setCSRData( targetNumRows, numColumns, targetIA, targetJA, targetValues );
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
    HArray<IndexType> sourceIA;

    HArray<IndexType> sourceJA;

    HArray<ValueType> sourceValues;

    other.buildCSRData( sourceIA, sourceJA, sourceValues );

    HArray<IndexType> targetIA;

    HArray<IndexType> targetJA;

    HArray<ValueType> targetValues;

    StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues, sourceIA, sourceJA, sourceValues,
            redistributor );

    const IndexType targetNumRows = targetIA.size() - 1;

    setCSRData( targetNumRows, numColumns, targetIA, targetJA, targetValues );
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
    HArray<IndexType> targetIA;

    HArray<IndexType> targetJA;

    HArray<ValueType> targetValues;

    StorageMethods<ValueType>::redistributeCSR( targetIA, targetJA, targetValues, other.getIA(), other.getJA(),
            other.getValues(), redistributor );

    const IndexType targetNumRows = targetIA.size() - 1;

    setCSRData( targetNumRows, numColumns, targetIA, targetJA, targetValues );
}

template<typename ValueType>
template<typename OtherValueType>
void MatrixStorage<ValueType>::setRawDenseData(
    const IndexType numRows,
    const IndexType numColumns,
    const OtherValueType values[] )
{
    // wrap all the data in a dense storage and make just an assign

    SCAI_LOG_INFO( logger, "set dense storage " << numRows << " x " << numColumns )

    HArrayRef<OtherValueType> denseData( numRows * numColumns, values );

    SCAI_LOG_INFO( logger, "use HArray array ref: " << denseData  )

    // move constructor avoids copy of the data

    DenseStorage<OtherValueType> denseStorage( numRows, numColumns, std::move( denseData ) );

    assign( denseStorage ); // will internally use the value epsilon

    SCAI_LOG_INFO( logger, *this << ": have set dense data " << numRows << " x " << numColumns )
}

/* ========================================================================= */
/*       File I/O                                                            */
/* ========================================================================= */

template<typename ValueType>
void MatrixStorage<ValueType>::writeToFile(
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType valuesType,
    const common::ScalarType indexType,
    const FileIO::FileMode mode ) const
{
    writeToFile( 1, 0, fileName, fileType, valuesType, indexType, mode );
}

template<typename ValueType>
void MatrixStorage<ValueType>::writeToFile(
    const PartitionId /* size */,
    const PartitionId /* rank */,
    const std::string& fileName,
    const std::string& fileType,
    const common::ScalarType dataType,
    const common::ScalarType indexType,
    const FileIO::FileMode mode ) const
{
    std::string suffix = fileType;

    if ( suffix == "" )
    {
        suffix = FileIO::getSuffix( fileName );
    }

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        if ( dataType != common::ScalarType::UNKNOWN )
        {
            fileIO->setDataType( dataType );
        }

        if ( indexType != common::ScalarType::UNKNOWN )
        {
            fileIO->setIndexType( indexType );
        }

        if ( mode != FileIO::DEFAULT_MODE )
        {
            fileIO->setMode( mode );
        }

        SCAI_LOG_INFO( logger, "write matrix storage to file, FileIO = " << *fileIO << ", storage = " << *this )

        // Note. SCAI_IO_TYPE_DATA allows that data is converted

        fileIO->writeStorage( *this, fileName );
    }
    else
    {
        COMMON_THROWEXCEPTION( "writeToFile " << fileName << ", unknown file type " << suffix )
    }
}

/*****************************************************************************/

template<typename ValueType>
void MatrixStorage<ValueType>::readFromFile( const std::string& fileName, const IndexType firstRow, IndexType nRows )
{
    SCAI_LOG_INFO( logger, "MatrixStorage<" << getValueType() << ">::readFromFile( " << fileName << ")" )
    SCAI_REGION( "Storage.readFromFile" )

    std::string suffix = FileIO::getSuffix( fileName );

    // Note: reading does not care about binary argument, just read as it is

    if ( FileIO::canCreate( suffix ) )
    {
        // okay, we can use FileIO class from factory

        std::unique_ptr<FileIO> fileIO( FileIO::create( suffix ) );

        // We do not set data type, take it from environment variable SCAI_IO_TYPE_DATA
        // fileIO->setDataType( common::ScalarType::INTERNAL );

        SCAI_LOG_INFO( logger, "Got from factory: " << *fileIO )

        fileIO->readStorage( *this, fileName, firstRow, nRows );
    }
    else
    {
        COMMON_THROWEXCEPTION( "readFromFile " << fileName << ", illegal suffix " << suffix )
    }

    check( "read matrix" );
}

/*****************************************************************************/

template<typename ValueType>
bool MatrixStorage<ValueType>::checkSymmetry() const
{
// check symmetry of matrix
    IndexType n = getNumRows();

    if ( n != getNumColumns() )
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

template<typename ValueType>
void MatrixStorage<ValueType>::fillCOO( 
    hmemo::HArray<IndexType> ia, 
    hmemo::HArray<IndexType> ja, 
    hmemo::HArray<ValueType> values,
    const common::BinaryOp op )
{
    // Default implementatinon convert this matrix to CSR, calls fillCOO for CSR storage, and convert bak

    SCAI_ASSERT_NE_ERROR( getFormat(), Format::CSR, "FATAL: CSR storage must override fillCOO" )

    IndexType numRows = getNumRows();
    IndexType numColumns = getNumColumns();

    HArray<IndexType> csrIA;
    HArray<IndexType> csrJA;
    HArray<ValueType> csrValues;

    buildCSRData( csrIA, csrJA, csrValues );

    CSRStorage<ValueType> csr( numRows, numColumns, std::move( csrIA ), std::move( csrJA ), std::move( csrValues ) );

    csr.fillCOO( std::move( ia ), std::move( ja ), std::move( values ), op );

    csr.splitUp( numRows, numColumns, csrIA, csrJA, csrValues );

    setCSRData( numRows, numColumns, csrIA, csrJA, csrValues );
}

/* ========================================================================= */
/*       Template Instantiations                                             */
/* ========================================================================= */

#define LAMA_MATRIXSTORAGE2_INST( ValueType, OtherValueType )                                          \
    template COMMON_DLL_IMPORTEXPORT void MatrixStorage<ValueType>::setRawDenseData<OtherValueType>(   \
            const IndexType, const IndexType, const OtherValueType* );                                 \

#define LAMA_MATRIXSTORAGE_INST( ValueType )                                                           \
    SCAI_COMMON_LOOP_LVL2( ValueType, LAMA_MATRIXSTORAGE2_INST, SCAI_NUMERIC_TYPES_HOST )

SCAI_COMMON_LOOP( LAMA_MATRIXSTORAGE_INST, SCAI_NUMERIC_TYPES_HOST )

#undef LAMA_MATRIXSTORAGE2_INST
#undef LAMA_MATRIXSTORAGE_INST


SCAI_COMMON_INST_CLASS( MatrixStorage, SCAI_NUMERIC_TYPES_HOST )

} /* end namespace lama */

} /* end namespace scai */
