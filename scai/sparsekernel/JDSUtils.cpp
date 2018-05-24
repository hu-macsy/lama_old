/**
 * @file JDSUtils.cpp
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
 * @brief Implementation and instantion of JDS utility methods.
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/JDSUtils.hpp>

#include <scai/sparsekernel/JDSKernelTrait.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>

namespace scai
{

using namespace hmemo;
using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( JDSUtils::logger, "CSRUtils" )

/* -------------------------------------------------------------------------- */

IndexType JDSUtils::getDiagonalPositions(
    HArray<IndexType>& diagonalPositions,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( numRows, jdsILG.size(), "illegally sized array jdsILG" )

    IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<JDSKernelTrait::getDiagonalPositions> kGetDiagonalPositions;

    // choose location where kernel routine is available

    ContextPtr loc = prefLoc;
    kGetDiagonalPositions.getSupportedContext( loc );

    ReadAccess<IndexType> rILG( jdsILG, loc );
    ReadAccess<IndexType> rDLG( jdsDLG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );

    WriteOnlyAccess<IndexType> wDiagonal( diagonalPositions, loc, numDiagonals );

    SCAI_CONTEXT_ACCESS( loc )

    IndexType numDiagonalsFound = 
        kGetDiagonalPositions[loc]( wDiagonal.get(), numDiagonals, numRows, 
                                    rILG.get(), rDLG.get(), rPerm.get(), rJA.get() );

    SCAI_LOG_INFO( logger, "getDiagonalPositions: " << numDiagonalsFound << " of " << numDiagonals << " available." )

    return numDiagonalsFound;
}

/* -------------------------------------------------------------------------- */

void JDSUtils::ilg2dlg(
    HArray<IndexType>& jdsDLG,
    const IndexType numDiagonals,
    const HArray<IndexType>& jdsILG,
    ContextPtr prefLoc )
{
    const IndexType numRows = jdsILG.size();

    // numDiagonals == jdsILG[0] 

    LAMAKernel<JDSKernelTrait::ilg2dlg> ilg2dlg;

    ContextPtr loc = prefLoc;

    ilg2dlg.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteOnlyAccess<IndexType> wDLG( jdsDLG, loc, numDiagonals );
    ReadAccess<IndexType> rILG( jdsILG, loc );
    ilg2dlg[loc]( wDLG.get(), numDiagonals, rILG.get(), numRows );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::getDiagonal(
    HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    ContextPtr prefLoc )
{
    HArray<IndexType> diagonalPositions;

    IndexType numDiagonalsFound = getDiagonalPositions( diagonalPositions, numRows, numColumns, 
                                                        jdsILG, jdsDLG, jdsPerm, jdsJA, prefLoc );

    // as we have the number of found diagonals we have not to check for any invalidIndex

    SCAI_ASSERT_EQ_ERROR( diagonalPositions.size(), numDiagonalsFound, 
                          "no diagonal property, some diagonal elements are missing" )

    HArrayUtils::gather( diagonal, jdsValues, diagonalPositions, common::BinaryOp::COPY, prefLoc );
}

/* -------------------------------------------------------------------------- */

#define JDSUTILS_SPECIFIER( ValueType )              \
                                                     \
    template void JDSUtils::getDiagonal(             \
            HArray<ValueType>&,                      \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \

SCAI_COMMON_LOOP( JDSUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef JDSUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
