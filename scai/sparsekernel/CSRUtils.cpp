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

IndexType CSRUtils::nonEmptyRows(        
    HArray<IndexType>& rowIndexes, 
    const HArray<IndexType>& csrIA,
    float threshold,
    ContextPtr )
{
    // currently only available on host

    ContextPtr loc = Context::getHostPtr();

    const IndexType numRows = csrIA.size() - 1;

    ReadAccess<IndexType> rIA( csrIA, loc );

    IndexType nonEmptyRows = OpenMPCSRUtils::countNonEmptyRowsByOffsets( rIA.get(), numRows );

    float usage = float( nonEmptyRows ) / float( numRows );

    if ( usage >= threshold )
    {
        SCAI_LOG_INFO( logger, "CSRStorage: do not build row indexes, usage = " << usage
                       << ", threshold = " << threshold )
    }
    else
    {
        SCAI_LOG_INFO( logger, "CSRStorage: build row indexes, #entries = " << nonEmptyRows )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, nonEmptyRows );
        OpenMPCSRUtils::setNonEmptyRowsByOffsets( wRowIndexes.get(), nonEmptyRows, rIA.get(), numRows );
    }

    return nonEmptyRows;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::compress( 
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const bool keepDiagonalValues,
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
        countNonZeros[loc]( wNewIA.get(), rIA.get(), rJA.get(), rValues.get(), numRows, eps, keepDiagonalValues );
    }

    newIA.resize( numRows );  //  reset size for scan1 operation

    // now compute the new offsets from the sizes, gives also new numValues

    IndexType newNumValues = HArrayUtils::scan1( newIA, prefContext );

    SCAI_LOG_ERROR( logger, "compress( keepDiagonals = " << keepDiagonalValues << ", eps = " << eps 
                     << " ) : "  << newNumValues << " non-diagonal zero elements, was " << numValues << " before" )

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
                           eps, keepDiagonalValues );
    }

    // now switch in place to the new data

    csrIA.swap( newIA );
    csrJA.swap( newJA );
    csrValues.swap( newValues );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::sort(
    HArray<IndexType>& ja,
    HArray<ValueType>& values,
    const HArray<IndexType>& ia,
    const IndexType numColumns,
    const bool diagonalFlag,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.CSR.sort" )

    SCAI_ASSERT_EQ_ERROR( values.size(), ja.size(), "serious size mismatch" )

    const IndexType numValues = values.size();   // #non-zero entries
    const IndexType numRows = ia.size() - 1;

    static LAMAKernel<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;

    ContextPtr loc = prefLoc;
    sortRowElements.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    WriteAccess<IndexType> wJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    sortRowElements[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numColumns, numValues, diagonalFlag );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType CSRUtils::setDiagonalFirst(
    HArray<IndexType>& ja,
    HArray<ValueType>& values,
    const HArray<IndexType>& ia,
    const IndexType numColumns,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.CSR.setDiagFirst" )

    const IndexType numRows = ia.size() - 1;

    static LAMAKernel<CSRKernelTrait::setDiagonalFirst<ValueType> > setDiagonal;

    ContextPtr loc = prefLoc;
    setDiagonal.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, "setDiagonalFirst on CSR data ( " << numRows << " x " << numColumns
                     << " ), called on " << *loc << ", preferred was " << *prefLoc )

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    WriteAccess<IndexType> wJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    IndexType numDiagonals = std::min( numRows, numColumns );

    return setDiagonal[loc]( wJA.get(), wValues.get(), numDiagonals, rIA.get() );
}


/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::convertCSR2CSC(
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

#define CSRUTILS_SPECIFIER( ValueType )              \
    template void CSRUtils::sort(                    \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const HArray<IndexType>&,                \
            const IndexType,                         \
            const bool,                              \
            ContextPtr );                            \
    template IndexType CSRUtils::setDiagonalFirst(   \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const HArray<IndexType>&,                \
            const IndexType,                         \
            ContextPtr );                            \
    template void CSRUtils::compress(                \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const bool,                              \
            const RealType<ValueType>,               \
            ContextPtr );                            \
    template void CSRUtils::convertCSR2CSC(          \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \

SCAI_COMMON_LOOP( CSRUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef CSRUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
