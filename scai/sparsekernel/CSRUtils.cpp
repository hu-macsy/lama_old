/**
 * @file CSRUtils.cpp
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
#include <scai/common/Constants.hpp>

#include <algorithm>

namespace scai
{

using namespace hmemo;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using tasking::SyncToken;

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

void CSRUtils::gatherSizes( 
    HArray<IndexType>& sizes,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& rowIndexes,
    ContextPtr prefLoc )
{
    LAMAKernel<CSRKernelTrait::gatherSizes> gatherSizes;

    ContextPtr loc = prefLoc;
    gatherSizes.getSupportedContext( loc );

    const IndexType nIndexes = rowIndexes.size();
    const IndexType numRows  = csrIA.size() - 1;

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rRows( rowIndexes, loc );
    WriteOnlyAccess<IndexType> wSizes( sizes, loc, nIndexes );
    gatherSizes[loc]( wSizes.get(), rIA.get(), numRows, rRows.get(), nIndexes );
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
    SCAI_ASSERT_EQ_DEBUG( rowJA.size(), rowValues.size(), "inconsistetne CSR data" )

    const IndexType numValues = rowJA.size();

    static LAMAKernel<CSRKernelTrait::convertCSR2CSC<ValueType> > convertCSR2CSC;
    ContextPtr loc = preferredLoc;
    convertCSR2CSC.getSupportedContext( loc );

    SCAI_LOG_DEBUG( logger,
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

template<typename ValueType>
RealType<ValueType> CSRUtils::absMaxDiffVal(
    const hmemo::HArray<IndexType>& aIA,
    const hmemo::HArray<IndexType>& aJA,
    const hmemo::HArray<ValueType>& aValues,
    const hmemo::HArray<IndexType>& bIA,
    const hmemo::HArray<IndexType>& bJA,
    const hmemo::HArray<ValueType>& bValues,
    const IndexType numRows,
    const IndexType numColumns,
    const bool isSorted,
    hmemo::ContextPtr prefLoc )
{
    if ( numRows == 0 || numColumns == 0 )
    {
        return static_cast<RealType<ValueType>>( 0 );
    }

    LAMAKernel<CSRKernelTrait::absMaxDiffVal<ValueType> > absMaxDiffVal;

    ContextPtr loc = prefLoc;
    absMaxDiffVal.getSupportedContext( loc );

    ReadAccess<IndexType> rCSRIA1( aIA, loc );
    ReadAccess<IndexType> rCSRJA1( aJA, loc );
    ReadAccess<ValueType> rCSRValues1( aValues, loc );
    ReadAccess<IndexType> rCSRIA2( bIA, loc );
    ReadAccess<IndexType> rCSRJA2( bJA, loc );
    ReadAccess<ValueType> rCSRValues2( bValues, loc );

    SCAI_CONTEXT_ACCESS( loc );

    auto maxVal = absMaxDiffVal[loc]( numRows, isSorted, 
                                      rCSRIA1.get(), rCSRJA1.get(), rCSRValues1.get(), 
                                      rCSRIA2.get(), rCSRJA2.get(), rCSRValues2.get() );
    return maxVal;
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

template<typename ValueType>
void CSRUtils::jacobiHalo(
    HArray<ValueType>& localSolution,
    const ValueType omega,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const HArray<IndexType>& rowIndexes,
    ContextPtr prefLoc )
{
    const IndexType numRows = localSolution.size();

    // not needed here: const IndexType numColumns = oldSolution.size();

    static LAMAKernel<CSRKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = prefLoc;

    jacobiHalo.getSupportedContext( loc );

    WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
    ReadAccess<ValueType> localDiagValues( localDiagonal, loc );
    ReadAccess<IndexType> haloIA( csrIA, loc );
    ReadAccess<IndexType> haloJA( csrJA, loc );
    ReadAccess<ValueType> haloValues( csrValues, loc );
    ReadAccess<ValueType> rOldSolution( oldHaloSolution, loc );

    const IndexType numNonEmptyRows = rowIndexes.size();

    if ( numNonEmptyRows != 0 )
    {
        SCAI_LOG_INFO( logger, "jacobiHalo optimized, #non-zero rows = " << numNonEmptyRows )

        ReadAccess<IndexType> rRowIndexes( rowIndexes, loc );
        SCAI_CONTEXT_ACCESS( loc )
        jacobiHalo[loc]( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                         rRowIndexes.get(), rOldSolution.get(), omega, numNonEmptyRows );
    }
    else
    {
        SCAI_LOG_INFO( logger, "jacobiHalo normal" )

        SCAI_CONTEXT_ACCESS( loc )
        jacobiHalo[loc]( wSolution.get(), localDiagValues.get(), haloIA.get(), haloJA.get(), haloValues.get(),
                         NULL, rOldSolution.get(), omega, numRows );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRUtils::gemv0(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // determine size of result as target of the linear mapping

    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), prefLoc );

        return NULL;   // already done
    }

    SCAI_REGION( "Sparse.CSR.gemv0" )

    ContextPtr loc = prefLoc;

    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );  

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL, 
                     numRows, numColumns, csrJA.size(), 
                     rIA.get(), rJA.get(), rValues.get(), op ); 
    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRUtils::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // if beta is 0, call the simpler routine, avoids accesses to y

    if ( beta == common::Constants::ZERO )
    {
        return gemv0( result, alpha, x, numRows, numColumns, csrIA, csrJA, csrValues, op, async, prefLoc );
    }

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        // result = beta * y, beta != 0

        HArrayUtils::compute( result, beta, common::BinaryOp::MULT, y, prefLoc );

        return NULL;
    }

    SCAI_REGION( "Sparse.CSR.gemv" )

    ContextPtr loc = prefLoc;

    static LAMAKernel<CSRKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, y.size() );  // okay if alias to y

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), 
                     numRows, numColumns, csrJA.size(), 
                     rIA.get(), rJA.get(), rValues.get(), op ); 
    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
tasking::SyncToken* CSRUtils::gemvSp(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const common::MatrixOp op,
    const HArray<IndexType>& nonZeroRowIndexes,
    bool async,
    ContextPtr prefLoc )
{
    if ( numRows == 0 || numColumns == 0 || alpha == 0 )
    {
        return NULL;
    }

    SCAI_REGION( "Sparse.CSR.gemvSp" )

    static LAMAKernel<CSRKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

    ContextPtr loc = prefLoc;

    sparseGEMV.getSupportedContext( loc );

    std::unique_ptr<tasking::SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    const IndexType numNonEmptyRows = nonZeroRowIndexes.size();

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );
    ReadAccess<IndexType> rIndexes( nonZeroRowIndexes, loc );

    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );

    sparseGEMV[loc]( wResult.get(),
                     alpha, rX.get(),
                     numNonEmptyRows, rIndexes.get(), rIA.get(), rJA.get(), rValues.get(), op );
    if ( async )
    {
        syncToken->pushRoutine( rIndexes.releaseDelayed() );
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRUtils::gemmSD0(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType k,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        HArrayUtils::setSameValue( result, nTarget * k, ValueType( 0 ), prefLoc );

        return NULL;   // already done
    }

    SCAI_LOG_INFO( logger, "gemmSD0: result[ " << numRows << " x " << k << " ] = " << alpha 
                          << " * csr [ " << numRows << " x " << numColumns << " ] * "
                          << " x [ " << numColumns << " x " << k << " ] " )

    ContextPtr loc = prefLoc;

    static LAMAKernel<CSRKernelTrait::gemmSD<ValueType> > gemmSD;

    gemmSD.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, nTarget * k );  

    gemmSD[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL, 
                 numRows, numColumns, k, 
                 rIA.get(), rJA.get(), rValues.get(), op ); 

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRUtils::gemmSD(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType k,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // if beta is 0, call the simpler routine, avoids accesses to y

    if ( beta == common::Constants::ZERO )
    {
        return gemmSD0( result, alpha, x, numRows, numColumns, k, csrIA, csrJA, csrValues, op, async, prefLoc );
    }

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        // result = beta * y, beta != 0

        HArrayUtils::compute( result, beta, common::BinaryOp::MULT, y, prefLoc );
        return NULL;
    }

    SCAI_LOG_INFO( logger, "gemmSD: result[ " << numRows << " x " << k << " ] = " << alpha 
                          << " * csr [ " << numRows << " x " << numColumns << " ] * "
                          << " x [ " << numColumns << " x " << k << " ] " 
                          << " + y [ " << numRows << " x " << k << " ] " )

    ContextPtr loc = prefLoc;

    static LAMAKernel<CSRKernelTrait::gemmSD<ValueType> > gemmSD;

    gemmSD.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, y.size() );  // okay if alias to y

    gemmSD[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), 
                 numRows, numColumns, k, 
                 rIA.get(), rJA.get(), rValues.get(), op ); 

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* CSRUtils::gemmDS(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType k,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "gemmDS: result[ " << k << " x " << numColumns << " ] = " << alpha 
                          << " * dense [ " << k << " x " << numRows << " ] * "
                          << " csr [ " << numRows << " x " << numColumns << " ] " 
                          << " + result [ " << k << " x " << numColumns << " ] " )

    ContextPtr loc = prefLoc;

    static LAMAKernel<CSRKernelTrait::gemmDS<ValueType> > gemmDS;

    gemmDS.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteAccess<ValueType> wResult( result ); 

    gemmDS[loc]( wResult.get(), alpha, rX.get(), beta, 
                 numRows, numColumns, k, 
                 rIA.get(), rJA.get(), rValues.get(), op ); 

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::reduce(
    HArray<ValueType>& result,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& csrIA,
    const hmemo::HArray<IndexType>& csrJA,
    const hmemo::HArray<ValueType>& csrValues,
    const IndexType dim,
    const common::BinaryOp reduceOp,
    const common::UnaryOp applyOp,
    hmemo::ContextPtr prefLoc )
{
    static LAMAKernel<CSRKernelTrait::reduce<ValueType> > reduce;

    ContextPtr loc = prefLoc;
    reduce.getSupportedContext( loc );

    IndexType resSize = dim == 0 ? numRows : numColumns;

    SCAI_ASSERT_EQ_ERROR( result.size(), resSize, "illegal size of result array" )

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );

    WriteAccess<ValueType> wResult( result, loc );

    reduce[loc]( wResult.get(),
                 rIA.get(), rJA.get(), rValues.get(), numRows, dim,
                 reduceOp, applyOp );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::solve(
    HArray<ValueType>& solution,
    const HArray<ValueType>& rhs,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& csrIA,
    const hmemo::HArray<IndexType>& csrJA,
    const hmemo::HArray<ValueType>& csrValues,
    bool isSymmetric,
    hmemo::ContextPtr prefLoc )
{
    if ( numRows == 0 || numColumns == 0 )
    {
        solution.clear();
        return;
    }

    LAMAKernel<CSRKernelTrait::decomposition<ValueType> > decomposition;

    ContextPtr loc = prefLoc;
    decomposition.getSupportedContext( loc );

    const IndexType numValues = csrJA.size();

    ReadAccess<IndexType> rCSRIA( csrIA, loc );
    ReadAccess<IndexType> rCSRJA( csrJA, loc );
    ReadAccess<ValueType> rCSRValues( csrValues, loc );
    ReadAccess<ValueType> rRHS( rhs, loc );
    WriteOnlyAccess<ValueType> wSol( solution, loc, numRows );
    SCAI_CONTEXT_ACCESS( loc );
    decomposition[loc]( wSol.get(), rCSRIA.get(), rCSRJA.get(), rCSRValues.get(),
                        rRHS.get(), numRows, numValues, isSymmetric );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::setRows(
    hmemo::HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& csrIA,
    const hmemo::HArray<IndexType>&,
    const hmemo::HArray<ValueType>& rowValues,
    const common::BinaryOp op,
    hmemo::ContextPtr prefLoc )
{
    if ( numRows == 0 || numColumns == 0 )
    {
        return;
    }

    static LAMAKernel<CSRKernelTrait::setRows<ValueType> > setRows;

    ContextPtr loc = prefLoc;

    setRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteAccess<ValueType> wValues( csrValues, loc );
    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<ValueType> rRows( rowValues, loc );

    setRows[loc]( wValues.get(), rIA.get(), numRows, rRows.get(), op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::setColumns(
    hmemo::HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& csrIA,
    const hmemo::HArray<IndexType>& csrJA,
    const hmemo::HArray<ValueType>& columnValues,
    const common::BinaryOp op,
    hmemo::ContextPtr prefLoc )
{
    if ( numRows == 0 || numColumns == 0 )
    {
        return;
    }

    static LAMAKernel<CSRKernelTrait::setColumns<ValueType> > setColumns;

    ContextPtr loc = prefLoc;

    setColumns.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteAccess<ValueType> wValues( csrValues, loc );
    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rCols( columnValues, loc );

    setColumns[loc]( wValues.get(), rIA.get(), rJA.get(), numRows, rCols.get(), op );
}

/* -------------------------------------------------------------------------- */

#define CSRUTILS_SPECIFIER( ValueType )                    \
                                                           \
    template void CSRUtils::sortRows(                      \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            ContextPtr );                                  \
                                                           \
    template IndexType CSRUtils::shiftDiagonalFirst(       \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::getDiagonal(                   \
            HArray<ValueType>&,                            \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const bool,                                    \
            ContextPtr );                                  \
                                                           \
    template bool CSRUtils::setDiagonalV(                  \
            HArray<ValueType>&,                            \
            const HArray<ValueType>&,                      \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const bool,                                    \
            ContextPtr );                                  \
                                                           \
    template bool CSRUtils::setDiagonal(                   \
            HArray<ValueType>&,                            \
            const ValueType,                               \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const bool,                                    \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::compress(                      \
            HArray<IndexType>&,                            \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const RealType<ValueType>,                     \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::convertCSR2CSC(                \
            HArray<IndexType>&,                            \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::matrixMultiply(                \
            HArray<IndexType>&,                            \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const ValueType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const IndexType,                               \
            const IndexType,                               \
            const IndexType,                               \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::matrixAdd(                     \
            HArray<IndexType>&,                            \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const ValueType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const ValueType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const IndexType,                               \
            const IndexType,                               \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::binaryOp(                      \
            HArray<IndexType>&,                            \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const IndexType,                               \
            const IndexType,                               \
            const common::BinaryOp,                        \
            ContextPtr );                                  \
                                                           \
    template RealType<ValueType> CSRUtils::absMaxDiffVal(  \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const IndexType,                               \
            const IndexType,                               \
            const bool,                                    \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::jacobi(                        \
            HArray<ValueType>&,                            \
            const ValueType,                               \
            const HArray<ValueType>&,                      \
            const HArray<ValueType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            ContextPtr );                                  \
                                                           \
    template void CSRUtils::jacobiHalo(                    \
            HArray<ValueType>&,                            \
            const ValueType,                               \
            const HArray<ValueType>&,                      \
            const HArray<ValueType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            const HArray<IndexType>&,                      \
            ContextPtr );                                  \
                                                           \
    template tasking::SyncToken* CSRUtils::gemv(           \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::MatrixOp,                            \
        const bool,                                        \
        ContextPtr );                                      \
                                                           \
    template tasking::SyncToken* CSRUtils::gemvSp(         \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::MatrixOp,                            \
        const HArray<IndexType>&,                          \
        const bool,                                        \
        ContextPtr );                                      \
                                                           \
    template tasking::SyncToken* CSRUtils::gemmSD(         \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const IndexType,                                   \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::MatrixOp,                            \
        const bool,                                        \
        ContextPtr );                                      \
                                                           \
    template tasking::SyncToken* CSRUtils::gemmDS(         \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const ValueType,                                   \
        const IndexType,                                   \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::MatrixOp,                            \
        const bool,                                        \
        ContextPtr );                                      \
                                                           \
    template void CSRUtils::reduce(                        \
        HArray<ValueType>&,                                \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const IndexType,                                   \
        const common::BinaryOp,                            \
        const common::UnaryOp,                             \
        ContextPtr );                                      \
                                                           \
    template void CSRUtils::setRows(                       \
        HArray<ValueType>&,                                \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::BinaryOp,                            \
        ContextPtr );                                      \
                                                           \
    template void CSRUtils::setColumns(                    \
        HArray<ValueType>&,                                \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::BinaryOp,                            \
        ContextPtr );                                      \
                                                           \
    template void CSRUtils::solve(                         \
        HArray<ValueType>&,                                \
        const HArray<ValueType>&,                          \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const bool,                                        \
        ContextPtr );                                      \

SCAI_COMMON_LOOP( CSRUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef CSRUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
