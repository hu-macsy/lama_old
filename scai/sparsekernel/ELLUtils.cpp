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

#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/SparseKernelTrait.hpp>
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
using utilskernel::SparseKernelTrait;

using sparsekernel::CSRUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( ELLUtils::logger, "ELLUtils" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::convertELL2CSR(
    hmemo::HArray<IndexType>& csrIA,
    hmemo::HArray<IndexType>& csrJA,
    hmemo::HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType,
    const hmemo::HArray<IndexType>& ellIA,
    const hmemo::HArray<IndexType>& ellJA,
    const hmemo::HArray<ValueType>& ellValues,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( ellIA.size(), numRows, "serious size mismatch" )

    if ( numRows == 0 )
    {
        csrIA = { 0 };
        csrJA.clear();
        csrValues.clear();
        return;
    }

    const IndexType numValuesPerRow = ellJA.size() / numRows;

    SCAI_ASSERT_EQ_DEBUG( numRows * numValuesPerRow, ellJA.size(), "serious size mismatch" )
    SCAI_ASSERT_EQ_DEBUG( ellValues.size(), ellJA.size(), "serious size mismatch" )

    const IndexType numValues = CSRUtils::sizes2offsets( csrIA, ellIA, prefLoc );

    // compute csrJA, csrValues

    static LAMAKernel<ELLKernelTrait::getCSRValues<ValueType> > getCSRValues;

    ContextPtr loc = prefLoc;

    getCSRValues.getSupportedContext( loc );

    ReadAccess<IndexType> rEllJA( ellJA, loc );
    ReadAccess<ValueType> rEllValues( ellValues, loc );
    ReadAccess<IndexType> rCsrIA( csrIA, loc );
    ReadAccess<IndexType> rEllIA( ellIA, loc );

    WriteOnlyAccess<IndexType> wCsrJA( csrJA, loc, numValues );
    WriteOnlyAccess<ValueType> wCsrValues( csrValues, loc, numValues );

    SCAI_CONTEXT_ACCESS( loc )

    getCSRValues[loc]( wCsrJA.get(), wCsrValues.get(), rCsrIA.get(), 
                       numRows, numValuesPerRow,
                       rEllIA.get(), rEllJA.get(), rEllValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::convertCSR2ELL(
    hmemo::HArray<IndexType>& ellIA,
    hmemo::HArray<IndexType>& ellJA,
    hmemo::HArray<ValueType>& ellValues,
    const IndexType numRows,
    const IndexType,
    const hmemo::HArray<IndexType>& csrIA,
    const hmemo::HArray<IndexType>& csrJA,
    const hmemo::HArray<ValueType>& csrValues,
    hmemo::ContextPtr prefLoc )
{
    CSRUtils::offsets2sizes( ellIA, csrIA, prefLoc );

    const IndexType numValuesPerRow = HArrayUtils::reduce( ellIA, common::BinaryOp::MAX, prefLoc );

    static LAMAKernel<ELLKernelTrait::setCSRValues<ValueType> > setCSRValues;

    ContextPtr loc = prefLoc;

    setCSRValues.getSupportedContext( loc );

    // now fill the array ellJA and ellValues

    ReadAccess<IndexType> rCsrIA( csrIA, loc );
    ReadAccess<IndexType> rCsrJA( csrJA, loc );
    ReadAccess<ValueType> rCsrValues( csrValues, loc );
    ReadAccess<IndexType> rEllIA( ellIA, loc );

    WriteOnlyAccess<IndexType> wEllJA( ellJA, loc, numRows * numValuesPerRow );
    WriteOnlyAccess<ValueType> wEllValues( ellValues, loc, numRows * numValuesPerRow );

    SCAI_LOG_DEBUG( logger, "convert CSR -> ELL, ellSize = " << numRows << " x " << numValuesPerRow )

    SCAI_CONTEXT_ACCESS( loc )

    setCSRValues[loc]( wEllJA.get(), wEllValues.get(), rEllIA.get(),
                       numRows, numValuesPerRow,
                       rCsrIA.get(), rCsrJA.get(), rCsrValues.get() );
}

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

        rowIndexes.clear();
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
    HArray<IndexType>& diagonalPositions,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "getDiagonalPositions for ellStrorage " << numRows << " x " << numColumns << ", size ja = " << ellJA.size() )

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
    HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    ContextPtr prefLoc )
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

template<typename ValueType>
void ELLUtils::setDiagonalV(
    HArray<ValueType>& ellValues,
    const HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    ContextPtr loc )
{
    HArray<IndexType> diagonalPositions;

    IndexType numDiagonalsFound = getDiagonalPositions( diagonalPositions, numRows, numColumns, ellIA, ellJA, loc );

    SCAI_ASSERT_EQ_ERROR( diagonal.size(), diagonalPositions.size(), "Illegally sized diagonal" )

    if ( diagonalPositions.size() == numDiagonalsFound )
    {
        bool unique = true;   // diagonal positons are unique for each row
        HArrayUtils::scatter( ellValues, diagonalPositions, unique, diagonal, common::BinaryOp::COPY, loc );
    }
    else
    {
        COMMON_THROWEXCEPTION( "cannot set diagonal, not all diagonal entries are available" )
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::setDiagonal(
    HArray<ValueType>& ellValues,
    const ValueType diagonalValue,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    ContextPtr loc )
{   
    HArray<IndexType> diagonalPositions;
    
    IndexType numDiagonalsFound = getDiagonalPositions( diagonalPositions, numRows, numColumns, ellIA, ellJA, loc );
    
    if ( diagonalPositions.size() == numDiagonalsFound )
    {   
        // we have no invalidIndex in diagonalPositions
        HArray<ValueType> diagonal( numDiagonalsFound, diagonalValue, loc );
        HArrayUtils::scatter( ellValues, diagonalPositions, true, diagonal, common::BinaryOp::COPY, loc );
    }
    else
    {
        COMMON_THROWEXCEPTION( "cannot set diagonal, not all diagonal entries are available" )
    }
}

/* -------------------------------------------------------------------------- */

IndexType ELLUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    ContextPtr prefLoc )
{   
    const IndexType numRows = ellIA.size();

    if ( numRows == 0 )
    {
        return invalidIndex;
    }

    const IndexType numValuesPerRow = ellJA.size() / numRows;

    // check row index to avoid out-of-range access, illegal j does not matter
    
    SCAI_ASSERT_VALID_INDEX_DEBUG( i, numRows, "row index out of range" )
    
    static LAMAKernel<ELLKernelTrait::getValuePos> getValuePos;
    
    ContextPtr loc = prefLoc;
    
    getValuePos.getSupportedContext( loc );
    
    SCAI_CONTEXT_ACCESS( loc )
    
    ReadAccess<IndexType> rIa( ellIA, loc );
    ReadAccess<IndexType> rJa( ellJA, loc );
    
    return getValuePos[loc]( i, j, numRows, numValuesPerRow, rIa.get(), rJa.get() );
}

/* -------------------------------------------------------------------------- */

void ELLUtils::getColumnPositions(
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>& positions,
    const hmemo::HArray<IndexType>& ellIA, 
    const HArray<IndexType>& ellJA,
    const IndexType j,
    const ContextPtr prefLoc )
{
    SCAI_REGION( "Storage.ELL.getSparseCol" )

    const IndexType numRows = ellIA.size();

    if ( numRows == 0 )
    {
        ia.clear();
        positions.clear();
        return;
    }

    const IndexType numValuesPerRow = ellJA.size() / numRows;    // numRows must not be 0

    static LAMAKernel<ELLKernelTrait::getColumnPositions> getColumnPositions;

    ContextPtr loc = prefLoc;

    getColumnPositions.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<IndexType> wRowIndexes( ia, loc, numRows );
    WriteOnlyAccess<IndexType> wValuePos( positions, loc, numRows );

    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<IndexType> rJA( ellJA, loc );

    IndexType cnt = getColumnPositions[loc]( wRowIndexes.get(), wValuePos.get(), j,
                                             rIA.get(), numRows, rJA.get(), numValuesPerRow );

    wRowIndexes.resize( cnt );
    wValuePos.resize( cnt );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::compress(
    HArray<IndexType>& ellIA,
    HArray<IndexType>& ellJA,
    HArray<ValueType>& ellValues,
    IndexType& numValuesPerRow, 
    const RealType<ValueType> eps,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.ELL.compress" )

    IndexType numRows = ellIA.size();

    SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ellJA.size(), "serious size mismatch for ELL arrays" )
    SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ellValues.size(), "serious size mismatch for ELL arrays" )

    HArray<IndexType> newIA;   // will contain size of each compressed row

    ContextPtr loc = prefLoc;

    static LAMAKernel<ELLKernelTrait::compressIA<ValueType> > compressIA;
    compressIA.getSupportedContext( loc );

    {
        ReadAccess<IndexType> rIA( ellIA, loc );
        ReadAccess<IndexType> rJA( ellJA, loc );
        ReadAccess<ValueType> rValues( ellValues, loc );

        WriteOnlyAccess<IndexType> wNewIA( newIA, loc, numRows );

        SCAI_CONTEXT_ACCESS( loc )

        compressIA[loc]( wNewIA.get(), rIA.get(), rJA.get(), rValues.get(), numRows, numValuesPerRow, eps );
    }

    IndexType nnzOld = HArrayUtils::reduce( ellIA, common::BinaryOp::ADD, prefLoc );
    IndexType nnzNew = HArrayUtils::reduce( newIA, common::BinaryOp::ADD, prefLoc );

    if ( nnzOld == nnzNew )
    {
        return;
    }

    IndexType newNumValuesPerRow = HArrayUtils::reduce( newIA, common::BinaryOp::MAX, prefLoc );

    static LAMAKernel<ELLKernelTrait::compressValues<ValueType> > compressValues;
    loc = prefLoc;
    compressValues.getSupportedContext( loc );

    if ( newNumValuesPerRow == numValuesPerRow )
    {
        // compress in place

        SCAI_CONTEXT_ACCESS( loc )

        ReadAccess<IndexType> rIA( ellIA, loc );
        WriteAccess<ValueType> wValues( ellValues, loc );
        WriteAccess<IndexType> wJA( ellJA, loc );
        compressValues[loc]( wJA.get(), wValues.get(), numValuesPerRow,
                             rIA.get(), wJA.get(), wValues.get(), numRows, numValuesPerRow, eps );
    }
    else
    {   // compress in new arrays

        HArray<ValueType> newValues;
        HArray<IndexType> newJA;

        IndexType newSize = numRows * newNumValuesPerRow;

        {
            SCAI_CONTEXT_ACCESS( loc )

            ReadAccess<IndexType> rOldIA( ellIA, loc );
            ReadAccess<IndexType> rOldJA( ellJA, loc );
            ReadAccess<ValueType> rOldValues( ellValues, loc );
            WriteOnlyAccess<ValueType> wNewValues( newValues, loc, newSize );
            WriteOnlyAccess<IndexType> wNewJA( newJA, loc, newSize );
            compressValues[loc]( wNewJA.get(), wNewValues.get(), newNumValuesPerRow,
                                 rOldIA.get(), rOldJA.get(), rOldValues.get(), numRows, numValuesPerRow, eps );
        }

        ellJA = std::move( newJA );
        ellValues = std::move( newValues );
        numValuesPerRow = newNumValuesPerRow;
    }

    ellIA = std::move( newIA );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* ELLUtils::gemv0(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // determine size of result as target of the linear mapping

    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    SCAI_ASSERT_EQ_DEBUG( ellValues.size(), ellJA.size(), "inconsistent sizes: ellJA - ellValues" )
    SCAI_ASSERT_EQ_DEBUG( numValuesPerRow * numRows, ellJA.size(), "inconsistent ellJA" )

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 || numValuesPerRow == 0 )
    {
        HArrayUtils::setSameValue( result, nTarget, ValueType( 0 ), prefLoc );

        return NULL;   // already done
    }

    ContextPtr loc = prefLoc;

    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<IndexType> rJA( ellJA, loc );
    ReadAccess<ValueType> rValues( ellValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL,
                     numRows, numColumns, numValuesPerRow,
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
SyncToken* ELLUtils::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    const common::MatrixOp op,
    const bool async,
    ContextPtr prefLoc )
{
    // if beta is 0, call the simpler routine, avoids accesses to y

    if ( beta == common::Constants::ZERO )
    {
        return gemv0( result, alpha, x, numRows, numColumns, numValuesPerRow, ellIA, ellJA, ellValues, op, async, prefLoc );
    }

    SCAI_ASSERT_EQ_DEBUG( numValuesPerRow * numRows, ellJA.size(), "inconsistent ellJA" )
    SCAI_ASSERT_EQ_DEBUG( ellValues.size(), ellJA.size(), "inconsistent ellJA - ellValues" )

    if ( alpha == common::Constants::ZERO  || numRows == 0 || numColumns == 0 || numValuesPerRow == 0 )
    {
        // result = beta * y, beta != 0

        SCAI_LOG_DEBUG( logger, "zero storage, set result = " << beta << " * y" )

        HArrayUtils::compute( result, beta, common::BinaryOp::MULT, y, prefLoc );

        return NULL;
    }

    SCAI_LOG_INFO( logger, "gemv: result = " << alpha << " * ELL " << op << " * x + " << beta << " * y " )

    SCAI_LOG_DEBUG( logger, "gemv: ELL is " << numRows << " x " << numColumns << ", #values/row = " << numValuesPerRow )

    ContextPtr loc = prefLoc;

    static LAMAKernel<ELLKernelTrait::normalGEMV<ValueType> > normalGEMV;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<IndexType> rJA( ellJA, loc );
    ReadAccess<ValueType> rValues( ellValues, loc );
    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );

    WriteOnlyAccess<ValueType> wResult( result, loc, y.size() );  // okay if alias to y

    normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(),
                     numRows, numColumns, numValuesPerRow,
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
tasking::SyncToken* ELLUtils::gemvSp(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const IndexType numRows,
    const IndexType numColumns,
    const IndexType numValuesPerRow,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    const common::MatrixOp op,
    const HArray<IndexType>& nonZeroRowIndexes,
    bool async,
    ContextPtr prefLoc )
{
    if ( numRows == 0 || numColumns == 0 || alpha == 0 || numValuesPerRow == 0 )
    {
        return NULL;
    }

    static LAMAKernel<ELLKernelTrait::sparseGEMV<ValueType> > sparseGEMV;

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

    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<IndexType> rJA( ellJA, loc );
    ReadAccess<ValueType> rValues( ellValues, loc );
    ReadAccess<IndexType> rIndexes( nonZeroRowIndexes, loc );

    ReadAccess<ValueType> rX( x, loc );
    WriteAccess<ValueType> wResult( result, loc );

    sparseGEMV[loc]( wResult.get(),
                     alpha, rX.get(),
                     numRows, numValuesPerRow, numNonEmptyRows, rIndexes.get(), 
                     rIA.get(), rJA.get(), rValues.get(), op );

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
tasking::SyncToken* ELLUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    const bool async,
    ContextPtr prefLoc )
{   
    SCAI_ASSERT_EQ_ERROR( rhs.size(), oldSolution.size(), "jacobi only for square matrices" )
    
    SCAI_ASSERT_EQ_ERROR( ellIA.size(), rhs.size(), "serious size mismatch for ELL arrays" ) 
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), ellValues.size(), "serious size mismatch for ELL arrays" )
    
    const IndexType numRows = rhs.size();

    if ( numRows == 0 )
    {
        return NULL;
    }

    const IndexType numColumns = oldSolution.size();
    const IndexType numValuesPerRow = ellJA.size() / numRows;
    
    static LAMAKernel<ELLKernelTrait::jacobi<ValueType>> jacobi;
    
    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );
    
    if ( &solution == &oldSolution )
    {   
        COMMON_THROWEXCEPTION( "alias of new/old solution is not allowed" )
    }
    
    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() )

    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<IndexType> rJA( ellJA, loc );
    ReadAccess<ValueType> rValues( ellValues, loc );
    
    ReadAccess<ValueType> rOld( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numColumns );
    
    SCAI_CONTEXT_ACCESS( loc );
    
    jacobi[ loc ]( wSolution.get(), numRows, numValuesPerRow, 
                   rIA.get(), rJA.get(), rValues.get(),
                   rOld.get(), rRhs.get(), omega );

    if ( async )
    {
        syncToken->pushRoutine( wSolution.releaseDelayed() );
        syncToken->pushRoutine( rRhs.releaseDelayed() );
        syncToken->pushRoutine( rOld.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::jacobiHalo(
    HArray<ValueType>& localSolution,
    const ValueType omega,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    const HArray<IndexType>& rowIndexes,
    ContextPtr prefLoc )
{
    const IndexType numRows = localSolution.size();
    const IndexType numValuesPerRow = ellValues.size() / numRows;

    SCAI_ASSERT_EQ_ERROR( numRows, ellIA.size(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ellJA.size(), "serious mismatch" )
    SCAI_ASSERT_EQ_ERROR( numRows * numValuesPerRow, ellValues.size(), "serious mismatch" )

    // not needed here: const IndexType numColumns = oldSolution.size();

    static LAMAKernel<ELLKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = prefLoc;

    jacobiHalo.getSupportedContext( loc );

    WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
    ReadAccess<ValueType> rDiagonal( localDiagonal, loc );
    ReadAccess<IndexType> haloIA( ellIA, loc );
    ReadAccess<IndexType> haloJA( ellJA, loc );
    ReadAccess<ValueType> haloValues( ellValues, loc );
    ReadAccess<ValueType> rOldSolution( oldHaloSolution, loc );

    const IndexType numNonEmptyRows = rowIndexes.size();

    if ( numNonEmptyRows != 0 )
    {
        SCAI_LOG_DEBUG( logger, "jacobiHalo optimized, #non-zero rows = " << numNonEmptyRows )

        ReadAccess<IndexType> rRowIndexes( rowIndexes, loc );
        SCAI_CONTEXT_ACCESS( loc )

        jacobiHalo[loc]( wSolution.get(), numRows, rDiagonal.get(), 
                         numValuesPerRow, haloIA.get(), haloJA.get(), haloValues.get(), 
                         rRowIndexes.get(), numNonEmptyRows, rOldSolution.get(), omega );
    }
    else
    {
        SCAI_LOG_DEBUG( logger, "jacobiHalo normal" )

        SCAI_CONTEXT_ACCESS( loc )

        jacobiHalo[loc]( wSolution.get(), numRows, rDiagonal.get(), 
                         numValuesPerRow, haloIA.get(), haloJA.get(), haloValues.get(), 
                         NULL, numRows, rOldSolution.get(), omega );
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void ELLUtils::setRows(
    hmemo::HArray<ValueType>& ellValues,
    const hmemo::HArray<IndexType>& ellIA,
    const hmemo::HArray<ValueType>& rowValues,
    const common::BinaryOp op,
    hmemo::ContextPtr prefLoc )
{
    const IndexType numRows = ellIA.size();

    if ( numRows == 0 || ellValues.size() == 0 )
    {
        return;
    }

    const IndexType numValuesPerRow = ellValues.size() / numRows;

    static LAMAKernel<ELLKernelTrait::setRows<ValueType> > setRows;

    ContextPtr loc = prefLoc;

    setRows.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    WriteAccess<ValueType> wValues( ellValues, loc );
    ReadAccess<IndexType> rIA( ellIA, loc );
    ReadAccess<ValueType> rRows( rowValues, loc );

    setRows[loc]( wValues.get(), numRows, numValuesPerRow, rIA.get(), rRows.get(), op );
}

/* -------------------------------------------------------------------------- */

#define ELLUTILS_SPECIFIER( ValueType )                    \
                                                           \
    template void ELLUtils::convertELL2CSR(                \
        HArray<IndexType>&,                                \
        HArray<IndexType>&,                                \
        HArray<ValueType>&,                                \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        ContextPtr );                                      \
                                                           \
    template void ELLUtils::convertCSR2ELL(                \
        HArray<IndexType>&,                                \
        HArray<IndexType>&,                                \
        HArray<ValueType>&,                                \
        const IndexType,                                   \
        const IndexType,                                   \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        ContextPtr );                                      \
                                                           \
    template void ELLUtils::getDiagonal(                   \
            HArray<ValueType>&,                            \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            const HArray<ValueType>&,                      \
            ContextPtr );                                  \
                                                           \
    template void ELLUtils::setDiagonalV(                  \
            HArray<ValueType>&,                            \
            const HArray<ValueType>&,                      \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            ContextPtr );                                  \
                                                           \
    template void ELLUtils::setDiagonal(                   \
            HArray<ValueType>&,                            \
            const ValueType,                               \
            const IndexType,                               \
            const IndexType,                               \
            const HArray<IndexType>&,                      \
            const HArray<IndexType>&,                      \
            ContextPtr );                                  \
                                                           \
    template void ELLUtils::compress(                      \
            HArray<IndexType>&,                            \
            HArray<IndexType>&,                            \
            HArray<ValueType>&,                            \
            IndexType&,                                    \
            const RealType<ValueType>,                     \
            ContextPtr );                                  \
                                                           \
    template tasking::SyncToken* ELLUtils::gemv(           \
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
    template tasking::SyncToken* ELLUtils::gemvSp(         \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const IndexType,                                   \
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
    template SyncToken* ELLUtils::jacobi(                  \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const HArray<ValueType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const bool,                                        \
        ContextPtr );                                      \
                                                           \
    template void ELLUtils::jacobiHalo(                    \
        HArray<ValueType>&,                                \
        const ValueType,                                   \
        const HArray<ValueType>&,                          \
        const HArray<ValueType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const HArray<IndexType>&,                          \
        ContextPtr );                                      \
                                                           \
    template void ELLUtils::setRows(                       \
        HArray<ValueType>&,                                \
        const HArray<IndexType>&,                          \
        const HArray<ValueType>&,                          \
        const common::BinaryOp,                            \
        ContextPtr );                                      \
      

SCAI_COMMON_LOOP( ELLUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef ELLUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
