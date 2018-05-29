/**
 * @file CSRUtils.cpp
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
 * @brief Implementation and instantion of CSR utility methods.
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>
#include <scai/sparsekernel/openmp/OpenMPCSRUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>

#include <algorithm>

namespace scai
{

using namespace hmemo;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( CSRUtils::logger, "CSRUtils" )

/* -------------------------------------------------------------------------- */

bool CSRUtils::validOffsets(
    const HArray<IndexType>& csrIA,
    const IndexType numValues,
    ContextPtr prefLoc )
{
    if ( csrIA.size() == 0 )
    {  
        return false;
    }

    LAMAKernel<CSRKernelTrait::validOffsets> validOffsets;

    ContextPtr loc = prefLoc;

    validOffsets.getSupportedContext( loc );

    IndexType numRows = csrIA.size() - 1;

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );

    return validOffsets[loc]( rIA.get(), numRows, numValues );
}

/* -------------------------------------------------------------------------- */

IndexType CSRUtils::sizes2offsets(
    HArray<IndexType>& offsets,
    const HArray<IndexType>& sizes,
    ContextPtr loc )
{
    const IndexType n = sizes.size();

    if ( &sizes != &offsets )
    {
        offsets.clear();
        offsets.reserve( loc, n +  1 );
        utilskernel::HArrayUtils::assign( offsets, sizes, loc );
    }

    return utilskernel::HArrayUtils::scan1( offsets, loc );
}

/* -------------------------------------------------------------------------- */

void CSRUtils::offsets2sizes (
    HArray<IndexType>& sizes,
    const HArray<IndexType>& offsets,
    ContextPtr )
{
    const IndexType n = offsets.size() - 1;

    if ( &sizes == &offsets )
    {
        auto wArray = hostWriteAccess( sizes );

        //   { 0, 5, 7, 11, 16  } -> { 5  2  4  5 }, no parallelism here

        for ( IndexType i = 0; i < n; i++ )
        {
            wArray[i] = wArray[i + 1] - wArray[i];
        }
 
        wArray.resize( n );
    }
    else
    { 
        auto rOffsets = hostReadAccess( offsets );
        auto wSizes = hostWriteOnlyAccess( sizes, n );

        for ( IndexType i = 0; i < n; i++ )
        {
            wSizes[i] = rOffsets[i + 1] - rOffsets[i];
        }
    }
}

/* -------------------------------------------------------------------------- */

IndexType CSRUtils::nonEmptyRows(        
    HArray<IndexType>& rowIndexes, 
    const HArray<IndexType>& csrIA,
    float threshold,
    ContextPtr prefLoc )
{
    // currently only available on host

    static LAMAKernel<CSRKernelTrait::nonEmptyRows> kNonEmptyRows;

    ContextPtr loc = prefLoc;

    kNonEmptyRows.getSupportedContext( loc );

    const IndexType numRows = csrIA.size() - 1;

    ReadAccess<IndexType> rIA( csrIA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType count = kNonEmptyRows[loc]( (IndexType*) NULL, rIA.get(), numRows );

    float usage = float( count ) / float( numRows );

    if ( usage >= threshold )
    {
        SCAI_LOG_INFO( logger, "CSRStorage: do not build row indexes, usage = " << usage
                       << ", threshold = " << threshold )

        rowIndexes.clear();
    }
    else
    {
        SCAI_LOG_INFO( logger, "CSRStorage: build row indexes, #entries = " << count )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, count );
        kNonEmptyRows[loc]( wRowIndexes.get(), rIA.get(), numRows );
    }

    return count;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::compress( 
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const RealType<ValueType> eps,
    ContextPtr prefContext )
{
    SCAI_REGION( "Sparse.CSR.compress" )

    static LAMAKernel<CSRKernelTrait::countNonZeros<ValueType> > countNonZeros;

    const IndexType numRows = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    HArray<IndexType> newIA;  // starts with sizes before it becomes offset array

    // compute the new sizes array

    {
        ContextPtr loc = prefContext;
        countNonZeros.getSupportedContext( loc );
        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<IndexType> rJA( csrJA, loc );
        ReadAccess<ValueType> rValues( csrValues, loc );
        WriteOnlyAccess<IndexType> wNewIA( newIA, loc, numRows + 1 );  // allocate already for offsets
        countNonZeros[loc]( wNewIA.get(), rIA.get(), rJA.get(), rValues.get(), numRows, eps );
    }

    newIA.resize( numRows );  //  reset size for scan1 operation

    // now compute the new offsets from the sizes, gives also new numValues

    IndexType newNumValues = HArrayUtils::scan1( newIA, prefContext );

    SCAI_LOG_INFO( logger, "compress( eps = " << eps << " ) : "  << newNumValues 
                           << " non-diagonal zero elements, was " << numValues << " before" )

    // ready if there are no new non-zero values

    if ( newNumValues == numValues )
    {
        return;
    }

    // All information is available how to fill the compressed data

    HArray<ValueType> newValues;
    HArray<IndexType> newJA;

    {
        static LAMAKernel<CSRKernelTrait::compress<ValueType> > compressData;
        ContextPtr loc = prefContext;
        compressData.getSupportedContext( loc );
        SCAI_CONTEXT_ACCESS( loc )
        ReadAccess<IndexType> rNewIA( newIA, loc );
        ReadAccess<IndexType> rIA( csrIA, loc );
        ReadAccess<IndexType> rJA( csrJA, loc );
        ReadAccess<ValueType> rValues( csrValues, loc );
        WriteOnlyAccess<IndexType> wNewJA( newJA, loc, newNumValues );
        WriteOnlyAccess<ValueType> wNewValues( newValues, loc, newNumValues );

        compressData[loc]( wNewJA.get(), wNewValues.get(), rNewIA.get(),
                           rIA.get(), rJA.get(), rValues.get(), numRows,
                           eps );
    }

    // now switch in place to the new data

    csrIA.swap( newIA );
    csrJA.swap( newJA );
    csrValues.swap( newValues );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::sortRows(
    HArray<IndexType>& ja,
    HArray<ValueType>& values,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.CSR.sort" )

    SCAI_ASSERT_EQ_DEBUG( values.size(), ja.size(), "serious size mismatch" )
    SCAI_ASSERT_EQ_DEBUG( numRows, ia.size() - 1, "serious size mismatch" )

    const IndexType numValues = values.size();   // #non-zero entries

    static LAMAKernel<CSRKernelTrait::sortRows<ValueType> > sortRows;

    ContextPtr loc = prefLoc;
    sortRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    WriteAccess<IndexType> wJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    sortRows[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numColumns, numValues );
}

/* -------------------------------------------------------------------------- */

bool CSRUtils::hasSortedRows(
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    hmemo::ContextPtr prefLoc )
{
    const IndexType numValues = ja.size();   // #non-zero entries

    static LAMAKernel<CSRKernelTrait::hasSortedRows> hasSortedRows;

    ContextPtr loc = prefLoc;
    hasSortedRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );

    bool is = hasSortedRows[loc]( rIA.get(), rJA.get(), numRows, numColumns, numValues );

    return is;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType CSRUtils::shiftDiagonalFirst(
    HArray<IndexType>& ja,
    HArray<ValueType>& values,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.CSR.shiftDiagFirst" )

    static LAMAKernel<CSRKernelTrait::shiftDiagonal<ValueType> > shiftDiagonal;

    ContextPtr loc = prefLoc;
    shiftDiagonal.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, "shiftDiagonalFirst on CSR data ( " << numRows << " x " << numColumns
                     << " ), called on " << *loc << ", preferred was " << *prefLoc )

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    WriteAccess<IndexType> wJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    return shiftDiagonal[loc]( wJA.get(), wValues.get(), numDiagonals, rIA.get() );
}


/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::convertCSR2CSC(
    HArray<IndexType>& colIA,
    HArray<IndexType>& colJA,
    HArray<ValueType>& colValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& rowIA,
    const HArray<IndexType>& rowJA,
    const HArray<ValueType>& rowValues,
    const ContextPtr preferredLoc )
{
    const IndexType numValues = rowJA.size();
    SCAI_ASSERT_EQUAL_DEBUG( rowJA.size(), rowValues.size() )
    static LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = preferredLoc;
    convertCSR2CSC.getSupportedContext( loc );
    SCAI_LOG_INFO( logger,
                   "MatrixStorage::CSR2CSC of matrix " << numRows << " x " << numColumns << ", #nnz = " << numValues << " on " << *loc )
    SCAI_REGION( "Sparse.CSR.2CSC" )
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

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::matrixMultiply(
    hmemo::HArray<IndexType>& cIA,
    hmemo::HArray<IndexType>& cJA,
    hmemo::HArray<ValueType>& cValues,
    const ValueType alpha,
    const hmemo::HArray<IndexType>& aIA,
    const hmemo::HArray<IndexType>& aJA,
    const hmemo::HArray<ValueType>& aValues,
    const hmemo::HArray<IndexType>& bIA,
    const hmemo::HArray<IndexType>& bJA,
    const hmemo::HArray<ValueType>& bValues,
    const IndexType m,
    const IndexType n,
    const IndexType k ,
    const hmemo::ContextPtr prefLoc )
{
    if ( m == 0 || n == 0 )
    {
        cIA = HArray<IndexType>( m + 1, ValueType( 0 ) );
        cJA.clear();
        cValues.clear();
        return;
    }

    static LAMAKernel<CSRKernelTrait::matrixMultiplySizes> matrixMultiplySizes;
    static LAMAKernel<CSRKernelTrait::matrixMultiply<ValueType> > matrixMultiply;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    matrixMultiply.getSupportedContext( loc, matrixMultiplySizes );

    ReadAccess<IndexType> rAIA( aIA, loc );
    ReadAccess<IndexType> rAJA( aJA, loc );
    ReadAccess<ValueType> rAValues( aValues, loc );
    ReadAccess<IndexType> rBIA( bIA, loc );
    ReadAccess<IndexType> rBJA( bJA, loc );
    ReadAccess<ValueType> rBValues( bValues, loc );

    WriteOnlyAccess<IndexType> wCIA( cIA, loc, m + 1 );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType nnz = matrixMultiplySizes[loc] ( wCIA.get(), m, n, k,
                                               rAIA.get(), rAJA.get(),
                                               rBIA.get(), rBJA.get() );

    WriteOnlyAccess<IndexType> wCJA( cJA, loc, nnz );
    WriteOnlyAccess<ValueType> wCValues( cValues, loc, nnz );

    matrixMultiply[loc]( wCIA.get(), wCJA.get(), wCValues.get(), m, n, k, alpha, 
                         rAIA.get(), rAJA.get(), rAValues.get(), rBIA.get(), rBJA.get(), rBValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::matrixAdd(
    hmemo::HArray<IndexType>& cIA,
    hmemo::HArray<IndexType>& cJA,
    hmemo::HArray<ValueType>& cValues,
    const ValueType alpha,
    const hmemo::HArray<IndexType>& aIA,
    const hmemo::HArray<IndexType>& aJA,
    const hmemo::HArray<ValueType>& aValues,
    const ValueType beta,
    const hmemo::HArray<IndexType>& bIA,
    const hmemo::HArray<IndexType>& bJA,
    const hmemo::HArray<ValueType>& bValues,
    const IndexType m,
    const IndexType n,
    const hmemo::ContextPtr prefLoc )
{
    if ( m == 0 || n == 0 )
    {
        cIA = HArray<IndexType>( m + 1, ValueType( 0 ) );
        cJA.clear();
        cValues.clear();
        return;
    }

    static LAMAKernel<CSRKernelTrait::matrixAddSizes> matrixAddSizes;
    static LAMAKernel<CSRKernelTrait::matrixAdd<ValueType> > matrixAdd;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    matrixAdd.getSupportedContext( loc, matrixAddSizes );

    ReadAccess<IndexType> rAIA( aIA, loc );
    ReadAccess<IndexType> rAJA( aJA, loc );
    ReadAccess<ValueType> rAValues( aValues, loc );
    ReadAccess<IndexType> rBIA( bIA, loc );
    ReadAccess<IndexType> rBJA( bJA, loc );
    ReadAccess<ValueType> rBValues( bValues, loc );

    WriteOnlyAccess<IndexType> wCIA( cIA, loc, m + 1 );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType nnz = matrixAddSizes[loc] ( wCIA.get(), m, n, 
                                          rAIA.get(), rAJA.get(),
                                          rBIA.get(), rBJA.get() );

    WriteOnlyAccess<IndexType> wCJA( cJA, loc, nnz );
    WriteOnlyAccess<ValueType> wCValues( cValues, loc, nnz );

    matrixAdd[loc]( wCJA.get(), wCValues.get(), wCIA.get(), 
                    m, n, 
                    alpha, rAIA.get(), rAJA.get(), rAValues.get(), 
                    beta,  rBIA.get(), rBJA.get(), rBValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::binaryOp(
    hmemo::HArray<IndexType>& cIA,
    hmemo::HArray<IndexType>& cJA,
    hmemo::HArray<ValueType>& cValues,
    const hmemo::HArray<IndexType>& aIA,
    const hmemo::HArray<IndexType>& aJA,
    const hmemo::HArray<ValueType>& aValues,
    const hmemo::HArray<IndexType>& bIA,
    const hmemo::HArray<IndexType>& bJA,
    const hmemo::HArray<ValueType>& bValues,
    const IndexType m,
    const IndexType n,
    const common::BinaryOp op,
    const hmemo::ContextPtr prefLoc )
{
    if ( m == 0 || n == 0 )
    {
        cIA = HArray<IndexType>( m + 1, ValueType( 0 ) );
        cJA.clear();
        cValues.clear();
        return;
    }

    SCAI_ASSERT_ERROR( hasSortedRows( m, n, aIA, aJA, prefLoc ), "binaryOp: input storage a not sorted" );
    SCAI_ASSERT_ERROR( hasSortedRows( m, n, bIA, bJA, prefLoc ), "binaryOp: input storage b not sorted" );

    static LAMAKernel<CSRKernelTrait::binaryOpSizes> binaryOpSizes;
    static LAMAKernel<CSRKernelTrait::binaryOp<ValueType> > binaryOp;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    binaryOp.getSupportedContext( loc, binaryOpSizes );

    ReadAccess<IndexType> rAIA( aIA, loc );
    ReadAccess<IndexType> rAJA( aJA, loc );
    ReadAccess<ValueType> rAValues( aValues, loc );

    ReadAccess<IndexType> rBIA( bIA, loc );
    ReadAccess<IndexType> rBJA( bJA, loc );
    ReadAccess<ValueType> rBValues( bValues, loc );

    WriteOnlyAccess<IndexType> wCIA( cIA, loc, m + 1 );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType nnz = binaryOpSizes[loc] ( wCIA.get(), m, n,
                                         rAIA.get(), rAJA.get(),
                                         rBIA.get(), rBJA.get() );

    WriteOnlyAccess<IndexType> wCJA( cJA, loc, nnz );
    WriteOnlyAccess<ValueType> wCValues( cValues, loc, nnz );

    binaryOp[loc]( wCJA.get(), wCValues.get(), wCIA.get(), 
                   m, n, 
                   rAIA.get(), rAJA.get(), rAValues.get(), 
                   rBIA.get(), rBJA.get(), rBValues.get(), op );
}

/* -------------------------------------------------------------------------- */

bool CSRUtils::hasDiagonalProperty(
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const bool isSorted,
    hmemo::ContextPtr prefLoc )
{
    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<CSRKernelTrait::hasDiagonalProperty> kHasDiagonalProperty;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    kHasDiagonalProperty.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );

    SCAI_CONTEXT_ACCESS( loc )

    return kHasDiagonalProperty[loc]( numDiagonals, rIA.get(), rJA.get(), isSorted );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::getDiagonal(
    hmemo::HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const hmemo::HArray<ValueType>& values,
    const bool isSorted,
    hmemo::ContextPtr prefLoc )
{
    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<CSRKernelTrait::getDiagonal<ValueType>> kGetDiagonal;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    kGetDiagonal.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    ReadAccess<ValueType> rValues( values, loc );

    WriteOnlyAccess<ValueType> wDiagonal( diagonal, loc, numDiagonals );

    SCAI_CONTEXT_ACCESS( loc )

    kGetDiagonal[loc]( wDiagonal.get(), numDiagonals, rIA.get(), rJA.get(), rValues.get(), isSorted );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
bool CSRUtils::setDiagonalV(
    hmemo::HArray<ValueType>& values,
    const hmemo::HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const bool isSorted,
    const hmemo::ContextPtr prefLoc )
{
    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    SCAI_ASSERT_EQ_ERROR( diagonal.size(), numDiagonals, "illegal size of diagonal to set for storage" )

    static LAMAKernel<CSRKernelTrait::setDiagonalV<ValueType> > setDiagonalV;

    // choose Context where all kernel routines are available
    
    ContextPtr loc = prefLoc;
    setDiagonalV.getSupportedContext( loc );
    
    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    ReadAccess<ValueType> rDiagonal( diagonal, loc );
    
    SCAI_CONTEXT_ACCESS( loc )

    bool okay = setDiagonalV[loc]( wValues.get(), rDiagonal.get(), numDiagonals, rIA.get(), rJA.get(), isSorted );

    return okay;
}

template<typename ValueType>
bool CSRUtils::setDiagonal(
    hmemo::HArray<ValueType>& values,
    const ValueType diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const bool isSorted,
    hmemo::ContextPtr prefLoc )
{
    IndexType numDiagonals = common::Math::min( numRows, numColumns );
    
    static LAMAKernel<CSRKernelTrait::setDiagonal<ValueType>> setDiagonal;

    // choose Context where all kernel routines are available
    
    ContextPtr loc = prefLoc;
    setDiagonal.getSupportedContext( loc );
    
    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    SCAI_CONTEXT_ACCESS( loc )

    bool okay = setDiagonal[loc]( wValues.get(), diagonal, numDiagonals, rIA.get(), rJA.get(), isSorted );

    return okay;
}

/* -------------------------------------------------------------------------- */

IndexType CSRUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    ContextPtr prefLoc )
{
    // check row index to avoid out-of-range access, illegal j does not matter

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, csrIA.size() - 1, "row index out of range" )

    static LAMAKernel<CSRKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = prefLoc;

    getValuePos.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( csrIA, loc );
    ReadAccess<IndexType> rJa( csrJA, loc );

    return getValuePos[loc]( i, j, rIa.get(), rJa.get() );
}

/* -------------------------------------------------------------------------- */

void CSRUtils::getColumnPositions(
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>& positions,
    const hmemo::HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const IndexType j,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.CSR.getColPos" )

    const IndexType numRows   = csrIA.size() - 1;
    const IndexType numValues = csrJA.size();

    static LAMAKernel<CSRKernelTrait::getColumnPositions> getColumnPositions;

    ContextPtr loc = prefLoc;

    getColumnPositions.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<IndexType> wRowIndexes( ia, loc, numRows );
    WriteOnlyAccess<IndexType> wValuePos( positions, loc, numRows );

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );

    IndexType cnt = getColumnPositions[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), numRows, rJA.get(), numValues );

    wRowIndexes.resize( cnt );
    wValuePos.resize( cnt );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( rhs.size(), oldSolution.size(), "jacobi only for square matrices" )

    SCAI_ASSERT_EQ_ERROR( csrIA.size() - 1, rhs.size(), "serious size mismatch for CSR arrays" )
    SCAI_ASSERT_EQ_ERROR( csrJA.size(), csrValues.size(), "serious size mismatch for CSR arrays" )

    const IndexType numRows = rhs.size();
    const IndexType numColumns = oldSolution.size();

    static LAMAKernel<CSRKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of new/old solution is not allowed" )
    }

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );

    ReadAccess<ValueType> rOld( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numColumns );

    SCAI_CONTEXT_ACCESS( loc );

    jacobi[loc]( wSolution.get(),
                 rIA.get(), rJA.get(), rValues.get(),
                 rOld.get(), rRhs.get(), omega, numRows );
}

/* -------------------------------------------------------------------------- */

#define CSRUTILS_SPECIFIER( ValueType )              \
                                                     \
    template void CSRUtils::sortRows(                \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            ContextPtr );                            \
                                                     \
    template IndexType CSRUtils::shiftDiagonalFirst( \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::getDiagonal(             \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const bool,                              \
            ContextPtr );                            \
                                                     \
    template bool CSRUtils::setDiagonalV(            \
            HArray<ValueType>&,                      \
            const HArray<ValueType>&,                \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const bool,                              \
            ContextPtr );                            \
                                                     \
    template bool CSRUtils::setDiagonal(             \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const bool,                              \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::compress(                \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const RealType<ValueType>,               \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::convertCSR2CSC(          \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::matrixMultiply(          \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const IndexType,                         \
            const IndexType,                         \
            const IndexType,                         \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::matrixAdd(               \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const ValueType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const IndexType,                         \
            const IndexType,                         \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::binaryOp(                \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            const IndexType,                         \
            const IndexType,                         \
            const common::BinaryOp,                  \
            ContextPtr );                            \
                                                     \
    template void CSRUtils::jacobi(                  \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const HArray<ValueType>&,                \
            const HArray<ValueType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \

SCAI_COMMON_LOOP( CSRUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef CSRUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
