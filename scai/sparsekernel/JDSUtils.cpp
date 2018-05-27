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
    SCAI_REGION( "Sparse.JDS.getDiagPos" )

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

IndexType JDSUtils::ilg2dlg(
    HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsILG,
    ContextPtr prefLoc )
{
    const IndexType numRows = jdsILG.size();

    if ( 0 == numRows )
    {
        jdsDLG.clear();
        return 0;
    }

    // numDiagonals, becomes size of array DLG, 

    IndexType numDiagonals = jdsILG[0];  // is maximal number of entries in row

    LAMAKernel<JDSKernelTrait::ilg2dlg> ilg2dlg;

    ContextPtr loc = prefLoc;

    ilg2dlg.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteOnlyAccess<IndexType> wDLG( jdsDLG, loc, numDiagonals );
    ReadAccess<IndexType> rILG( jdsILG, loc );

    return ilg2dlg[loc]( wDLG.get(), numDiagonals, rILG.get(), numRows );
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

IndexType JDSUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    ContextPtr prefLoc )
{
    const IndexType numRows = jdsILG.size();

    // check row index to avoid out-of-range access, illegal j does not matter

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, numRows, "row index out of range" )

    static LAMAKernel<JDSKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = prefLoc;

    getValuePos.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rJa( jdsJA, loc );

    IndexType pos = getValuePos[loc]( i, j, numRows, rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

    return pos;
}

/* -------------------------------------------------------------------------- */

void JDSUtils::getColumnPositions(
    HArray<IndexType>& ia,
    HArray<IndexType>& positions,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const IndexType j,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.JDS.getColPos" )

    const IndexType numRows = jdsILG.size();

    static LAMAKernel<JDSKernelTrait::getColumnPositions> getColumnPositions;

    ContextPtr loc = prefLoc;

    getColumnPositions.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<IndexType> wRowIndexes( ia, loc, numRows );
    WriteOnlyAccess<IndexType> wValuePos( positions, loc, numRows );

    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rJa( jdsJA, loc );

    IndexType cnt = getColumnPositions[loc]( wRowIndexes.get(), wValuePos.get(), j, numRows,
                                             rIlg.get(), rDlg.get(), rPerm.get(), rJa.get() );

    wRowIndexes.resize( cnt );
    wValuePos.resize( cnt );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.JDS.jacobi" )

    SCAI_ASSERT_EQ_ERROR( rhs.size(), oldSolution.size(), "jacobi only for square matrices" )

    SCAI_ASSERT_EQ_ERROR( jdsILG.size(), rhs.size(), "serious size mismatch for JDS arrays" )
    SCAI_ASSERT_EQ_ERROR( jdsJA.size(), jdsValues.size(), "serious size mismatch for JDS arrays" )

    const IndexType numRows = rhs.size();
    const IndexType numDiagonals = jdsDLG.size();

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of new/old solution is not allowed" )
    }

    static LAMAKernel<JDSKernelTrait::jacobi<ValueType>> jacobi;

    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );

    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    ReadAccess<ValueType> rValues( jdsValues, loc );

    ReadAccess<ValueType> rOld( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numRows );

    SCAI_CONTEXT_ACCESS( loc );

    jacobi[loc]( wSolution.get(), numRows, rPerm.get(), 
                 rIlg.get(), numDiagonals, rDlg.get(), rJA.get(),
                 rValues.get(), rOld.get(), rRhs.get(), omega );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::jacobiHalo(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& diagonal,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.JDS.jacobiHalo" )

    SCAI_ASSERT_EQ_ERROR( jdsJA.size(), jdsValues.size(), "serious size mismatch for JDS arrays" )

    SCAI_ASSERT_EQ_ERROR( diagonal.size(), jdsILG.size(), "serious size mismatch for JDS arrays/diagonal" )
    SCAI_ASSERT_EQ_ERROR( jdsPerm.size(), jdsILG.size(), "serious size mismatch for JDS arrays" )

    const IndexType numRows = jdsILG.size();
    const IndexType numDiagonals = jdsDLG.size();

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of new/old solution is not allowed" )
    }

    static LAMAKernel<JDSKernelTrait::jacobiHalo<ValueType>> jacobiHalo;

    ContextPtr loc = prefLoc;
    jacobiHalo.getSupportedContext( loc );

    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    ReadAccess<ValueType> rValues( jdsValues, loc );

    ReadAccess<ValueType> rOld( oldSolution, loc );
    ReadAccess<ValueType> rDiagonal( diagonal, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numRows );

    SCAI_CONTEXT_ACCESS( loc );

    jacobiHalo[loc]( wSolution.get(), numRows, rDiagonal.get(), numDiagonals,
                     rPerm.get(), rIlg.get(), rDlg.get(), rJA.get(),
                     rValues.get(), rOld.get(), omega );
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
                                                     \
    template void JDSUtils::jacobi(                  \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const HArray<ValueType>&,                \
            const HArray<ValueType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \
                                                     \
    template void JDSUtils::jacobiHalo(              \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const HArray<ValueType>&,                \
            const HArray<ValueType>&,                \
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
