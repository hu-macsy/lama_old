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
    const IndexType numValuesPerRow = ellJA.size() / numRows;

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
void ELLUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const HArray<IndexType>& ellIA,
    const HArray<IndexType>& ellJA,
    const HArray<ValueType>& ellValues,
    ContextPtr prefLoc )
{   
    SCAI_ASSERT_EQ_ERROR( rhs.size(), oldSolution.size(), "jacobi only for square matrices" )
    
    SCAI_ASSERT_EQ_ERROR( ellIA.size(), rhs.size(), "serious size mismatch for ELL arrays" ) 
    SCAI_ASSERT_EQ_ERROR( ellJA.size(), ellValues.size(), "serious size mismatch for CSR arrays" )
    
    const IndexType numRows = rhs.size();
    const IndexType numColumns = oldSolution.size();
    const IndexType numValuesPerRow = ellJA.size() / numRows;
    
    static LAMAKernel<ELLKernelTrait::jacobi<ValueType>> jacobi;
    
    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );
    
    if ( &solution == &oldSolution )
    {   
        COMMON_THROWEXCEPTION( "alias of new/old solution is not allowed" )
    }
    
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
    template void ELLUtils::setDiagonalV(            \
            HArray<ValueType>&,                      \
            const HArray<ValueType>&,                \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            ContextPtr );                            \
    template void ELLUtils::setDiagonal(             \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const IndexType,                         \
            const IndexType,                         \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            ContextPtr );                            \
    template void ELLUtils::compress(                \
            HArray<IndexType>&,                      \
            HArray<IndexType>&,                      \
            HArray<ValueType>&,                      \
            IndexType&,                              \
            const RealType<ValueType>,               \
            ContextPtr );                            \
    template void ELLUtils::jacobi(                  \
            HArray<ValueType>&,                      \
            const ValueType,                         \
            const HArray<ValueType>&,                \
            const HArray<ValueType>&,                \
            const HArray<IndexType>&,                \
            const HArray<IndexType>&,                \
            const HArray<ValueType>&,                \
            ContextPtr );                            \

SCAI_COMMON_LOOP( ELLUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef ELLUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
