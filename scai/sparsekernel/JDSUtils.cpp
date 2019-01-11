/**
 * @file JDSUtils.cpp
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
 * @brief Implementation and instantion of JDS utility methods.
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/JDSUtils.hpp>

#include <scai/sparsekernel/JDSKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>
#include <scai/common/Constants.hpp>

namespace scai
{

using namespace hmemo;

using tasking::SyncToken;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( JDSUtils::logger, "JDSUtils" ) 

/* -------------------------------------------------------------------------- */

void JDSUtils::buildRowSizes(
    HArray<IndexType>& rowSizes,
    const HArray<IndexType>& jdsIlg,
    const HArray<IndexType>& jdsPerm,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( jdsIlg.size(), jdsPerm.size(), "serious mismatch" )

    const IndexType numRows = jdsIlg.size();

    rowSizes.clear();
    rowSizes.reserve( prefLoc, numRows + 1 );  // avoid reallocation when building offsets

    rowSizes.resize( numRows );             // no initialization required

    HArrayUtils::setSameValue<IndexType>( rowSizes, numRows, 0, prefLoc );

    HArrayUtils::scatter( rowSizes, jdsPerm, true, jdsIlg, common::BinaryOp::COPY, prefLoc );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::convertJDS2CSR(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& jdsIlg,
    const HArray<IndexType>& jdsDlg,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( jdsIlg.size(), numRows, "JDS ilg array illege size" )

    buildRowSizes( csrIA, jdsIlg, jdsPerm, prefLoc );

    IndexType numValues = HArrayUtils::scan1( csrIA );   // build  the offset aray

    SCAI_LOG_INFO( logger, "convertJDS2CSR of storage " << numRows << " x " << numColumns << ", nnz = " << numValues )

    SCAI_ASSERT_EQ_ERROR( numValues, jdsValues.size(), "row sizes do not sum up to number of nnz entries" )

    SCAI_LOG_DEBUG( logger, "buildCSR from JDS with " << numValues << " values" )

    // compute the inverse permutation so that we find original row in JDS data

    HArray<IndexType> jdsInvPerm; 
    HArrayUtils::inversePerm( jdsInvPerm, jdsPerm, prefLoc);

    static LAMAKernel<JDSKernelTrait::getCSRValues<ValueType>> getCSRValues;

    ContextPtr loc = prefLoc;

    getCSRValues.getSupportedContext( loc );

    ReadAccess<IndexType> rJdsILG( jdsIlg, loc );
    ReadAccess<IndexType> rJdsDLG( jdsDlg, loc );
    ReadAccess<IndexType> rJdsInversePerm( jdsInvPerm, loc );
    ReadAccess<IndexType> rJdsJA( jdsJA, loc );
    ReadAccess<ValueType> rJdsValues( jdsValues, loc );

    ReadAccess<IndexType> rCsrIA( csrIA, loc );
    WriteOnlyAccess<IndexType> wCsrJA( csrJA, loc, numValues ); 
    WriteOnlyAccess<ValueType> wCsrValues( csrValues, loc, numValues );

    SCAI_CONTEXT_ACCESS( loc )

    // now we can convert JDS to CSR by JDS kernel routine

    getCSRValues[loc]( wCsrJA.get(), wCsrValues.get(), rCsrIA.get(), 
                       numRows, rJdsInversePerm.get(), rJdsILG.get(),
                       rJdsDLG.get(), rJdsJA.get(), rJdsValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::convertCSR2JDS(
    HArray<IndexType>& jdsIlg,
    HArray<IndexType>& jdsDlg,
    HArray<IndexType>& jdsPerm,
    HArray<IndexType>& jdsJA,
    HArray<ValueType>& jdsValues,
    const IndexType numRows,
    const IndexType,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    ContextPtr prefLoc )
{
    // fill up the array ilg and perm

    CSRUtils::offsets2sizes( jdsIlg, csrIA, prefLoc );
    HArrayUtils::setOrder( jdsPerm, numRows, prefLoc );

    // sort ilg in descending order, keep permutation

    HArrayUtils::sort( &jdsPerm, &jdsIlg, jdsIlg, false, prefLoc );

    // build dlg array by ilg

    IndexType numValues = ilg2dlg( jdsDlg, jdsIlg, prefLoc );  // dlg is now correctly defined

    SCAI_ASSERT_EQ_ERROR( numValues, csrJA.size(), "sum of row sizes does not match size of ja" )
    SCAI_ASSERT_EQ_ERROR( numValues, csrValues.size(), "sum of row sizes does not match size of values." )

    IndexType numDiagonals = jdsDlg.size();

    static LAMAKernel<JDSKernelTrait::setCSRValues<ValueType>> setCSRValues;

    ContextPtr loc = prefLoc;
    setCSRValues.getSupportedContext( loc );

    ReadAccess<IndexType> rCsrIA( csrIA, loc );
    ReadAccess<IndexType> rCsrJA( csrJA, loc );
    ReadAccess<ValueType> rCsrValues( csrValues, loc );

    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rIlg( jdsIlg, loc );
    ReadAccess<IndexType> rDlg( jdsDlg, loc );

    WriteOnlyAccess<ValueType> wValues( jdsValues, loc, numValues );
    WriteOnlyAccess<IndexType> wJa( jdsJA, loc, numValues );

    SCAI_CONTEXT_ACCESS( loc )

    setCSRValues[loc]( wJa.get(), wValues.get(), numRows, rPerm.get(), rIlg.get(), numDiagonals, rDlg.get(),
                       rCsrIA.get(), rCsrJA.get(), rCsrValues.get() ); 
}

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

void JDSUtils::getRowPositions(
    HArray<IndexType>& positions,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const IndexType i,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.JDS.getColPos" )

    SCAI_ASSERT_EQ_DEBUG( jdsILG.size(), jdsPerm.size(), "jds ilg and perm have inconsistent sizes" )

    const IndexType numRows = jdsILG.size();

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, numRows, "row index out of range" )

    const IndexType numDiagonals = jdsDLG.size();  // maximal entries in one row

    static LAMAKernel<JDSKernelTrait::getRowPositions> getRowPositions;

    ContextPtr loc = prefLoc;

    getRowPositions.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<IndexType> wPositions( positions, loc, numDiagonals );

    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );

    IndexType cnt = getRowPositions[loc]( wPositions.get(), i, numRows,
                                          rIlg.get(), rDlg.get(), rPerm.get() );

    SCAI_ASSERT_LE_DEBUG( cnt, numDiagonals, "serious error: inconsistent JDS" )

    wPositions.resize( cnt );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::setRow(
    HArray<ValueType>& jdsValues,
    const IndexType i,
    const HArray<ValueType>& row,
    const HArray<IndexType>& jdsIlg,
    const HArray<IndexType>& jdsDlg,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const common::BinaryOp op,
    ContextPtr prefLoc )
{
    const IndexType numColumns = row.size();
    const IndexType numRows = jdsIlg.size();

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, numRows, "row index out of range" )

    SCAI_LOG_INFO( logger, "setRowImpl( i = " << i << " )" )

    static LAMAKernel<JDSKernelTrait::setRow<ValueType> > setRow;

    ContextPtr loc = prefLoc;

    setRow.getSupportedContext( loc );

    ReadAccess<IndexType> rIlg( jdsIlg, loc );
    ReadAccess<IndexType> rDlg( jdsDlg, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    WriteAccess<ValueType> wValues( jdsValues, loc );
    ReadAccess<ValueType> rRow( row, loc );
    SCAI_CONTEXT_ACCESS( loc )
    setRow[loc]( wValues.get(), i, numColumns, numRows, 
                 rPerm.get(), rIlg.get(), rDlg.get(), rJA.get(), rRow.get(), op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::getRow(
    HArray<ValueType>& row,
    const IndexType numColumns,
    const IndexType i,
    const HArray<IndexType>& jdsIlg,
    const HArray<IndexType>& jdsDlg,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    ContextPtr prefLoc )
{   
    const IndexType numRows = jdsIlg.size();

    SCAI_ASSERT_VALID_INDEX_DEBUG( i, numRows, "row index out of range" )
    
    SCAI_LOG_INFO( logger, "setRowImpl( i = " << i << " )" )
    
    static LAMAKernel<JDSKernelTrait::getRow<ValueType> > getRow;
    
    ContextPtr loc = prefLoc;

    getRow.getSupportedContext( loc );
    
    ReadAccess<IndexType> rIlg( jdsIlg, loc );
    ReadAccess<IndexType> rDlg( jdsDlg, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    ReadAccess<ValueType> rValues( jdsValues, loc );

    WriteOnlyAccess<ValueType> wRow( row, loc, numColumns );

    SCAI_CONTEXT_ACCESS( loc ) 

    getRow[loc]( wRow.get(), i, numColumns, numRows, rPerm.get(), rIlg.get(), rDlg.get(), rJA.get(), rValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* JDSUtils::gemv0(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& jdsIlg,
    const HArray<IndexType>& jdsDlg,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), prefLoc );

        return NULL;   // already done
    }

    const IndexType numDiagonals = jdsDlg.size();

    ContextPtr loc = prefLoc;

    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rDLG( jdsDlg, loc );
    ReadAccess<IndexType> rILG( jdsIlg, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    ReadAccess<ValueType> rValues( jdsValues, loc );

    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );  // okay if alias to y

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL,
                     numRows, numColumns, rPerm.get(), rILG.get(),
                     numDiagonals, rDLG.get(), rJA.get(), rValues.get(), op );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rPerm.releaseDelayed() );
        syncToken->pushRoutine( rDLG.releaseDelayed() );
        syncToken->pushRoutine( rILG.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* JDSUtils::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& jdsIlg,
    const HArray<IndexType>& jdsDlg,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // if beta is 0, call the simpler routine, avoids accesses to y

    if ( beta == common::Constants::ZERO )
    {
        return gemv0( result, alpha, x, numRows, numColumns, 
                      jdsIlg, jdsDlg, jdsPerm, jdsJA, jdsValues, op, async, prefLoc );
    }

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 )
    {
        // result = beta * y, beta != 0

        HArrayUtils::compute( result, beta, common::BinaryOp::MULT, y, prefLoc );

        return NULL;
    }

    ContextPtr loc = prefLoc;

    static LAMAKernel<JDSKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    const IndexType numDiagonals = jdsDlg.size();

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rDLG( jdsDlg, loc );
    ReadAccess<IndexType> rILG( jdsIlg, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    ReadAccess<ValueType> rValues( jdsValues, loc );

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, y.size() );  // okay if alias to y

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(),
                     numRows, numColumns, rPerm.get(), rILG.get(),
                     numDiagonals, rDLG.get(), rJA.get(), rValues.get(), op );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rPerm.releaseDelayed() );
        syncToken->pushRoutine( rDLG.releaseDelayed() );
        syncToken->pushRoutine( rILG.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* JDSUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& jdsValues,
    bool async,
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

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

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

    if ( async )
    {
        syncToken->pushRoutine( wSolution.releaseDelayed() );
        syncToken->pushRoutine( rRhs.releaseDelayed() );
        syncToken->pushRoutine( rOld.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rIlg.releaseDelayed() );
        syncToken->pushRoutine( rDlg.releaseDelayed() );
        syncToken->pushRoutine( rPerm.releaseDelayed() );
    }

    return syncToken.release();
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

template<typename ValueType>
void JDSUtils::setRows(
    HArray<ValueType>& jdsValues,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>&,
    const HArray<ValueType>& rowValues,
    const common::BinaryOp op,
    ContextPtr prefLoc )
{
    const IndexType numRows = jdsILG.size();

    if ( numRows == 0 || jdsValues.size() == 0 )
    {
        return;
    }

    static LAMAKernel<JDSKernelTrait::setRows<ValueType> > setRows;

    ContextPtr loc = prefLoc;

    setRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteAccess<ValueType> wValues( jdsValues, loc );
    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<ValueType> rRows( rowValues, loc );

    setRows[loc]( wValues.get(), numRows, rPerm.get(), rIlg.get(), rDlg.get(), rRows.get(), op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void JDSUtils::setColumns(
    HArray<ValueType>& jdsValues,
    const HArray<IndexType>& jdsILG,
    const HArray<IndexType>& jdsDLG,
    const HArray<IndexType>& jdsPerm,
    const HArray<IndexType>& jdsJA,
    const HArray<ValueType>& columnValues,
    const common::BinaryOp op,
    ContextPtr prefLoc )
{
    const IndexType numRows = jdsILG.size();

    if ( numRows == 0 || jdsValues.size() == 0 )
    {
        return;
    }

    static LAMAKernel<JDSKernelTrait::setColumns<ValueType> > setColumns;

    ContextPtr loc = prefLoc;

    setColumns.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteAccess<ValueType> wValues( jdsValues, loc );
    ReadAccess<IndexType> rIlg( jdsILG, loc );
    ReadAccess<IndexType> rDlg( jdsDLG, loc );
    ReadAccess<IndexType> rPerm( jdsPerm, loc );
    ReadAccess<IndexType> rJA( jdsJA, loc );
    ReadAccess<ValueType> rColumns( columnValues, loc );

    setColumns[loc]( wValues.get(), numRows, rPerm.get(), rIlg.get(), rDlg.get(), rJA.get(), rColumns.get(), op );
}

/* -------------------------------------------------------------------------- */

#define JDSUTILS_SPECIFIER( ValueType )              \
                                                     \
    template void JDSUtils::convertJDS2CSR(          \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::convertCSR2JDS(          \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::getDiagonal(             \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::getRow(                  \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::setRow(                  \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const common::BinaryOp,                      \
        ContextPtr );                                \
                                                     \
    template tasking::SyncToken* JDSUtils::gemv(     \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        const common::MatrixOp,                      \
        const bool,                                  \
        ContextPtr );                                \
                                                     \
    template SyncToken* JDSUtils::jacobi(            \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const HArray<ValueType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        const bool,                                  \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::jacobiHalo(              \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const HArray<ValueType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::setRows(                 \
        HArray<ValueType>&,                          \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        const common::BinaryOp,                      \
        ContextPtr );                                \
                                                     \
    template void JDSUtils::setColumns(              \
        HArray<ValueType>&,                          \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        const common::BinaryOp,                      \
        ContextPtr );                                \


SCAI_COMMON_LOOP( JDSUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef JDSUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
