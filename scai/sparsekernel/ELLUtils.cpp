/**
 * @file ELLUtils.cpp
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
 * @brief Implementation and instantion of ELL utility methods.
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/ELLUtils.hpp>

#include <scai/sparsekernel/ELLKernelTrait.hpp>

#include <scai/utilskernel/SparseKernelTrait.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>

namespace scai
{

using namespace hmemo;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using utilskernel::SparseKernelTrait;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( ELLUtils::logger, "CSRUtils" )

/* -------------------------------------------------------------------------- */

IndexType ELLUtils::nonEmptyRows(        
    HArray<IndexType>& rowIndexes, 
    const HArray<IndexType>& ellIA,
    float threshold,
    ContextPtr prefLoc )
{
    const IndexType numRows = ellIA.size();

    if ( numRows == 0 )
    {
        rowIndexes.clear();
        return 0;
    }

    static LAMAKernel<SparseKernelTrait::countNonZeros<IndexType> > countNonZeros;
    static LAMAKernel<SparseKernelTrait::compress<IndexType, IndexType> > compress;

    // choose location where both routines are available

    ContextPtr loc = prefLoc;
    countNonZeros.getSupportedContext( loc, compress );

    ReadAccess<IndexType> rIA( ellIA, loc );

    SCAI_CONTEXT_ACCESS( loc )

    // count the number of non-zero rows to have a good value for allocation of rowIndexes

    const IndexType ZERO = 0;   // sparse storage uses always the real 0
    const IndexType EPS  = 0;   // no tolerances used here

    IndexType nonEmptyRows = countNonZeros[loc]( rIA.get(), numRows, ZERO, EPS );

    float usage = float( nonEmptyRows ) / float( numRows );

    if ( usage >= threshold )
    {
        SCAI_LOG_INFO( logger, "ELLStorage: do not build row indexes, usage = " << usage
                       << ", threshold = " << threshold )
    }
    else
    {
        SCAI_LOG_INFO( logger, "ELLStorage: build row indexes, #entries = " << nonEmptyRows )

        WriteOnlyAccess<IndexType> wRowIndexes( rowIndexes, loc, nonEmptyRows );

        IndexType cnt = compress[loc]( NULL, wRowIndexes.get(), rIA.get(), numRows, ZERO, EPS );

        SCAI_ASSERT_EQ_ERROR( cnt, nonEmptyRows, "serious mismatch" );
    }

    return nonEmptyRows;
}

/* -------------------------------------------------------------------------- */

IndexType ELLUtils::getDiagonalPositions(
    hmemo::HArray<IndexType>& diagonalPositions,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ellIA,
    const hmemo::HArray<IndexType>& ellJA,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( numRows, ellIA.size(), "illegally sized array ellIA" )

    // Deal with case numRows == 0 just here as we divide later by numRows

    if ( numRows == 0 )
    {
        diagonalPositions.clear();
        return 0;
    }

    IndexType numDiagonals = common::Math::min( numRows, numColumns );
    IndexType numValuesPerRow = ellJA.size() / numRows;

    SCAI_ASSERT_EQ_DEBUG( numRows * numValuesPerRow, ellJA.size(), "illegal size" )

    static LAMAKernel<ELLKernelTrait::getDiagonalPositions> kGetDiagonalPositions;

    // choose location where kernel routine is available

    ContextPtr loc = prefLoc;
    kGetDiagonalPositions.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<IndexType> rJA( ellJA, loc );

    WriteOnlyAccess<IndexType> wDiagonal( diagonalPositions, loc, numDiagonals );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType numDiagonalsFound = 
        kGetDiagonalPositions[loc]( wDiagonal.get(), numDiagonals, numRows, numValuesPerRow, rIA.get(), rJA.get() );

    SCAI_LOG_INFO( logger, "getDiagonalPositions: " << numDiagonalsFound << " of " << numDiagonals << " available." )

    return numDiagonalsFound;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::getDiagonal(
    hmemo::HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ellIA,
    const hmemo::HArray<IndexType>& ellJA,
    const hmemo::HArray<ValueType>& ellValues,
    hmemo::ContextPtr prefLoc )
{
    HArray<IndexType> diagonalPositions;

    IndexType numDiagonalsFound = getDiagonalPositions( diagonalPositions, numRows, numColumns, ellIA, ellJA, prefLoc );

    // as we have the number of found diagonals we have not to check for any invalidIndex

    SCAI_ASSERT_EQ_ERROR( diagonalPositions.size(), numDiagonalsFound, 
                          "no diagonal property, some diagonal elements are missing" )

    HArrayUtils::gather( diagonal, ellValues, diagonalPositions, common::BinaryOp::COPY, prefLoc );

    // ToDo: getDiagonal might also be defined if not all diagonal entries are available
}

/* -------------------------------------------------------------------------- */

#define ELLUTILS_SPECIFIER( ValueType )              \
                                                     \
    template void ELLUtils::getDiagonal(             \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \

SCAI_COMMON_LOOP( ELLUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef ELLUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
