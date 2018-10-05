/**
 * @file COOUtils.cpp
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
 * @brief Utility functions for COO data
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/COOUtils.hpp>

#include <scai/sparsekernel/COOKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

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
using tasking::SyncToken;

namespace sparsekernel
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( COOUtils::logger, "COOUtils" )

/* -------------------------------------------------------------------------- */

void COOUtils::convertCOO2CSR(
    HArray<IndexType>& csrIA,
    const HArray<IndexType>& cooIA,
    const IndexType numRows,
    ContextPtr prefLoc )
{
    const IndexType numValues = cooIA.size();

    static LAMAKernel<COOKernelTrait::ia2offsets> ia2offsets;

    ContextPtr loc = prefLoc;

    ia2offsets.getSupportedContext( loc );

    WriteOnlyAccess<IndexType> wOffsets( csrIA, loc, numRows + 1 );
    ReadAccess<IndexType> rIA( cooIA, loc );
   
    SCAI_CONTEXT_ACCESS( loc )

    SCAI_LOG_INFO( logger, "ia2offsets, build csrIA[" << numRows << "] from cooIA[" << numValues << "] @ " << *loc )

    ia2offsets[loc]( wOffsets.get(), numRows, rIA.get(), numValues );
}

/* -------------------------------------------------------------------------- */

void COOUtils::convertCSR2COO(
    HArray<IndexType>& cooIA,
    const HArray<IndexType>& csrIA,
    const IndexType numValues,
    ContextPtr prefLoc )
{
    static LAMAKernel<COOKernelTrait::offsets2ia> offsets2ia;

    const IndexType numRows = csrIA.size() - 1;

    ContextPtr loc = prefLoc;

    offsets2ia.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    WriteOnlyAccess<IndexType> wIA( cooIA, loc, numValues );
    ReadAccess<IndexType> rOffsets( csrIA, loc );

    offsets2ia[loc]( wIA.get(), numValues, rOffsets.get(), numRows );
}

/* -------------------------------------------------------------------------- */

bool COOUtils::isSorted(
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    ContextPtr )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "inconsitent size of coo arrays" )

    const IndexType numValues = ia.size();

    if ( numValues < 2 )
    {
        return true;   // Important: numValues - 1 might be illegal value
    }

    bool sorted = true;

    auto rIA = hostReadAccess( ia );
    auto rJA = hostReadAccess( ja );

    #pragma omp parallel for

    for ( IndexType k = 0; k < numValues - 1; ++k )
    {
        if ( rIA[k] > rIA[k+1] )
        {
            sorted = false;
        }
        else if ( rIA[k] == rIA[k+1] && rJA[k] >= rJA[k + 1] )
        {
            // entries at same position are also not allowed
            sorted = false;
        }
    }
    return sorted;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::sort(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>& values,
    ContextPtr )
{
    SCAI_REGION( "Sparse.COO.sort" )

    using namespace utilskernel;

    const IndexType nnz = values.size();

    SCAI_ASSERT_EQ_ERROR( nnz, ia.size(), "illegal size for ia of COO" )
    SCAI_ASSERT_EQ_ERROR( nnz, ja.size(), "illegal size for ja of COO" )

    ContextPtr host = Context::getHostPtr();

    HArray<IndexType> perm;
    HArrayUtils::setOrder( perm, nnz, host );

    struct cmp
    {
        cmp( const HArray<IndexType>& ia, const HArray<IndexType>& ja )
        {
            pI = hostReadAccess( ia ).begin();
            pJ = hostReadAccess( ja ).begin();
        }

        bool operator()( IndexType p1, IndexType p2 )
        {
            if ( pI[p1] < pI[p2] )
            {
                return true;
            }
            else if ( pI[p1] > pI[p2] )
            {
                return false;
            }
            else
            {
                return pJ[p1] < pJ[p2];
            }
        }
 
        const IndexType *pI;
        const IndexType *pJ;
    };

    {
        cmp myCmp( ia, ja );
        auto wPerm = hostWriteAccess( perm );
        std::stable_sort( wPerm.begin(), wPerm.end() , myCmp );
    }

    HArray<IndexType> iaOld( std::move( ia ) );
    HArray<IndexType> jaOld( std::move( ja ) );
    HArray<ValueType> valuesOld( std::move( values ) );

    HArrayUtils::gather( ia, iaOld, perm, common::BinaryOp::COPY, host );
    HArrayUtils::gather( ja, jaOld, perm, common::BinaryOp::COPY, host );
    HArrayUtils::gather( values, valuesOld, perm, common::BinaryOp::COPY, host );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::unique(
    HArray<IndexType>& ia,
    HArray<IndexType>& ja,
    HArray<ValueType>& values,
    common::BinaryOp op,
    ContextPtr )

{
    SCAI_REGION( "Sparse.COO.unique" )

    const IndexType nnz = values.size();

    SCAI_ASSERT_EQ_ERROR( ia.size(), nnz, "serious mismatch for COO data" )
    SCAI_ASSERT_EQ_ERROR( ja.size(), nnz, "serious mismatch for COO data" )

    IndexType lastI = invalidIndex;
    IndexType lastJ = invalidIndex;

    IndexType nnzU  = 0;   // running index, ends with nnz for unique values

    {
        auto wIA = hostWriteAccess( ia );
        auto wJA = hostWriteAccess( ja );
        auto wValues = hostWriteAccess( values );

        for ( IndexType k = 0; k < nnz; ++k )
        {
            IndexType i = wIA[k];
            IndexType j = wJA[k];

            if ( i == lastI && j == lastJ )
            {
                wValues[nnzU - 1] = common::applyBinary( wValues[nnzU - 1], op, wValues[k] );
            }
            else
            {
                wValues[nnzU] = wValues[k];
                wIA[nnzU] = wIA[k];
                wJA[nnzU] = wJA[k];
                nnzU++;
            }

            lastI = i;
            lastJ = j;
        }
    }

    ia.resize( nnzU );
    ja.resize( nnzU );
    values.resize( nnzU );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::normalize(
    HArray<IndexType>& cooIA,
    HArray<IndexType>& cooJA,
    HArray<ValueType>& cooValues,
    common::BinaryOp op,
    ContextPtr prefLoc )
{
    if ( COOUtils::isSorted( cooIA, cooJA, prefLoc ) )
    {
        return;
    }

    sort( cooIA, cooJA, cooValues, prefLoc );
    unique( cooIA, cooJA, cooValues, op, prefLoc );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::getDiagonal(
    HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    const HArray<ValueType>& values,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( ia.size(), values.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = ia.size();
    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<COOKernelTrait::getDiagonal<ValueType>> kGetDiagonal;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    kGetDiagonal.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    ReadAccess<ValueType> rValues( values, loc );

    WriteOnlyAccess<ValueType> wDiagonal( diagonal, loc, numDiagonals );

    SCAI_CONTEXT_ACCESS( loc )

    kGetDiagonal[loc]( wDiagonal.get(), numDiagonals, rIA.get(), rJA.get(), rValues.get(), numValues );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::setDiagonalV(
    HArray<ValueType>& values,
    const HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( ia.size(), values.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = ia.size();
    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<COOKernelTrait::setDiagonalV<ValueType>> kSetDiagonalV;

    // choose context where kernel routines is available

    ContextPtr loc = prefLoc;
    kSetDiagonalV.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    ReadAccess<ValueType> rDiagonal( diagonal, loc );

    SCAI_CONTEXT_ACCESS( loc )

    kSetDiagonalV[loc]( wValues.get(), rDiagonal.get(), numDiagonals, rIA.get(), rJA.get(), numValues );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::setDiagonal(
    HArray<ValueType>& values,
    const ValueType diagonalValue,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( ia.size(), values.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = ia.size();
    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<COOKernelTrait::setDiagonal<ValueType>> kSetDiagonal;

    // choose context where kernel routines is available

    ContextPtr loc = prefLoc;
    kSetDiagonal.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    SCAI_CONTEXT_ACCESS( loc )

    kSetDiagonal[loc]( wValues.get(), diagonalValue, numDiagonals, rIA.get(), rJA.get(), numValues );
}

/* -------------------------------------------------------------------------- */

bool COOUtils::hasDiagonal(
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& ia,
    const HArray<IndexType>& ja,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = ia.size();
    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    if ( numDiagonals == 0 )
    {
        return true;
    }

    if ( numValues < numDiagonals )
    {
        return false;
    }

    static LAMAKernel<COOKernelTrait::hasDiagonalProperty> kHasDiagonalProperty;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    kHasDiagonalProperty.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );

    SCAI_CONTEXT_ACCESS( loc )

    bool hasIt = kHasDiagonalProperty[loc]( numDiagonals, rIA.get(), rJA.get(), numValues );

    return hasIt;
}

/* -------------------------------------------------------------------------- */

IndexType COOUtils::getValuePos( 
    const IndexType i, 
    const IndexType j, 
    const HArray<IndexType>& ia, 
    const HArray<IndexType>& ja,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( ia.size(), ja.size(), "inconsisent sizes of coo arrays." )

    IndexType numValues = ia.size();

    static LAMAKernel<COOKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = prefLoc;

    getValuePos.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIa( ia, loc );
    ReadAccess<IndexType> rJa( ja, loc );

    IndexType pos = getValuePos[loc]( i, j, rIa.get(), rJa.get(), numValues );

    return pos;
}

/* -------------------------------------------------------------------------- */

void COOUtils::getColumnPositions( 
    HArray<IndexType>& positions,
    const HArray<IndexType>& cooJA,
    const IndexType j,
    ContextPtr prefLoc )
{
    const IndexType numValues = cooJA.size();

    static LAMAKernel<COOKernelTrait::getColumn> getColumn;

    ContextPtr loc = prefLoc;

    getColumn.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rJA( cooJA, loc );

    // first call to determine the number of entries

    IndexType nColumnEntries = getColumn[loc]( NULL, rJA.get(), numValues, j );

    WriteOnlyAccess<IndexType> wPos( positions, loc, nColumnEntries );
    
    // second call to get the positions of the entries

    getColumn[loc]( wPos.get(), rJA.get(), numValues, j );
}

/* -------------------------------------------------------------------------- */

void COOUtils::getRowPositions(
    IndexType& offset,
    IndexType& n,
    const HArray<IndexType>& cooIA,
    const IndexType i,
    ContextPtr prefLoc )
{
    const IndexType numValues = cooIA.size();

    static LAMAKernel<COOKernelTrait::getRow> getRow;

    ContextPtr loc = prefLoc;

    getRow.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( cooIA, loc );

    n = getRow[loc]( offset, rIA.get(), numValues, i );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::setRows(
    HArray<ValueType>& cooValues,
    const IndexType numRows, 
    const IndexType,
    const HArray<IndexType>& cooIA,
    const HArray<IndexType>&,
    const HArray<ValueType>& rowValues,
    const common::BinaryOp op,         
    ContextPtr prefLoc )
{
    const IndexType numValues = cooValues.size();

    // Note: the column indexes of the entries do not matter here at all

    SCAI_ASSERT_EQ_ERROR( numValues, cooIA.size(), "inconsistent COO storage data" )
    SCAI_ASSERT_EQ_ERROR( numRows, rowValues.size(), "illegal size for array with update values for each row" )

    static LAMAKernel<COOKernelTrait::setRows<ValueType>> setRows;

    ContextPtr loc = prefLoc;

    setRows.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rRowValues( rowValues, loc );
    WriteAccess<ValueType> wValues( cooValues, loc );  // update
    ReadAccess<IndexType> rIA( cooIA, loc );
    setRows[loc]( wValues.get(), rIA.get(), rRowValues.get(), numValues, op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::setColumns(
    HArray<ValueType>& cooValues,
    const IndexType,
    const IndexType numColumns,
    const HArray<IndexType>&,
    const HArray<IndexType>& cooJA,
    const HArray<ValueType>& columnValues,
    const common::BinaryOp op,         
    ContextPtr prefLoc )
{
    const IndexType numValues = cooValues.size();

    // Note: the row indexes of the entries do not matter here at all

    SCAI_ASSERT_EQ_ERROR( numValues, cooJA.size(), "inconsistent COO storage data" )
    SCAI_ASSERT_EQ_ERROR( numColumns, columnValues.size(), "illegal size for array with update values for each column" )

    static LAMAKernel<COOKernelTrait::setColumns<ValueType>> setColumns;
 
    ContextPtr loc = prefLoc;

    setColumns.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<ValueType> rColumnValues( columnValues, loc );
    WriteAccess<ValueType> wValues( cooValues, loc );  // update
    ReadAccess<IndexType> rJA( cooJA, loc );
    setColumns[loc]( wValues.get(), rJA.get(), rColumnValues.get(), numValues, op );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::binaryOp(
    HArray<IndexType>& cIA,
    HArray<IndexType>& cJA,
    HArray<ValueType>& cValues,
    const HArray<IndexType>& aIA,
    const HArray<IndexType>& aJA,
    const HArray<ValueType>& aValues,
    const HArray<IndexType>& bIA,
    const HArray<IndexType>& bJA,
    const HArray<ValueType>& bValues,
    const IndexType m,
    const IndexType n,
    const common::BinaryOp op,
    const ContextPtr prefLoc )
{
    if ( m == 0 || n == 0 )
    {
        cIA.clear();
        cJA.clear();
        cValues.clear();
        return;
    }

    // convert COO ia arrays to offsets array and use CSR utility for binaryOp
    // Allows also for better OpenMP/CUDA parallelization

    HArray<IndexType> aOffsets;
    HArray<IndexType> bOffsets;
    HArray<IndexType> cOffsets;

    convertCOO2CSR( aOffsets, aIA, m, prefLoc );
    convertCOO2CSR( bOffsets, bIA, m, prefLoc );

    CSRUtils::binaryOp( cOffsets, cJA, cValues, 
                        aOffsets, aJA, aValues,
                        bOffsets, bJA, bValues,
                        m, n, op, prefLoc );

    convertCSR2COO( cIA, cOffsets, cJA.size(), prefLoc );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::jacobi(
    HArray<ValueType>& solution,
    const ValueType omega,
    const HArray<ValueType>& oldSolution,
    const HArray<ValueType>& rhs,
    const HArray<IndexType>& cooIA,
    const HArray<IndexType>& cooJA,
    const HArray<ValueType>& cooValues,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( rhs.size(), oldSolution.size(), "jacobi only for square matrices" )
   
    const IndexType numRows = rhs.size();
    const IndexType numColumns = oldSolution.size();

    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooJA.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooValues.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = cooIA.size();

    static LAMAKernel<COOKernelTrait::jacobi<ValueType> > jacobi;

    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );

    if ( &solution == &oldSolution )
    {
        COMMON_THROWEXCEPTION( "alias of new/old solution is not allowed" )
    }

    ReadAccess<IndexType> rIA( cooIA, loc );
    ReadAccess<IndexType> rJA( cooJA, loc );
    ReadAccess<ValueType> rValues( cooValues, loc );

    ReadAccess<ValueType> rOld( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numColumns );

    SCAI_CONTEXT_ACCESS( loc );

    jacobi[loc]( wSolution.get(),
                 numValues, rIA.get(), rJA.get(), rValues.get(),
                 rOld.get(), rRhs.get(), omega, numRows );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::jacobiHalo(
    HArray<ValueType>& localSolution,
    const ValueType omega,
    const HArray<ValueType>& localDiagonal,
    const HArray<ValueType>& oldHaloSolution,
    const HArray<IndexType>& cooIA,
    const HArray<IndexType>& cooJA,
    const HArray<ValueType>& cooValues,
    ContextPtr prefLoc )
{
    const IndexType numRows = localSolution.size();
    const IndexType numValues = cooIA.size();

    SCAI_ASSERT_EQ_DEBUG( numValues, cooJA.size(), "inconsistent COO data" )
    SCAI_ASSERT_EQ_DEBUG( numValues, cooValues.size(), "inconsistent COO data" )
    SCAI_ASSERT_EQ_DEBUG( numRows, localDiagonal.size(), "diaogonal elements for each row required" )

    if ( numRows == 0 || numValues == 0 )
    {
        return;   // avoids overhead of accesses + kernel calls
    }

    // not needed here: const IndexType numColumns = oldSolution.size();

    static LAMAKernel<COOKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    ContextPtr loc = prefLoc;

    jacobiHalo.getSupportedContext( loc );

    WriteAccess<ValueType> wSolution( localSolution, loc ); // will be updated
    ReadAccess<ValueType> localDiagValues( localDiagonal, loc );
    ReadAccess<IndexType> haloIA( cooIA, loc );
    ReadAccess<IndexType> haloJA( cooJA, loc );
    ReadAccess<ValueType> haloValues( cooValues, loc );
    ReadAccess<ValueType> rOldSolution( oldHaloSolution, loc );

    SCAI_CONTEXT_ACCESS( loc )

    jacobiHalo[loc]( wSolution.get(), 
                     numValues, haloIA.get(), haloJA.get(), haloValues.get(),
                     localDiagValues.get(), rOldSolution.get(), omega, numRows );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* COOUtils::gemv(
    HArray<ValueType>& result,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& cooIA,
    const HArray<IndexType>& cooJA,
    const HArray<ValueType>& cooValues,
    const common::MatrixOp op,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    bool async,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.COO.gemv" )

    // SCAI_ASSERT_EQ_ERROR( x.size(), numColumns, "illegally sized array x for gemv" )
    // SCAI_ASSERT_EQ_ERROR( y.size(), numRows, "illegally sized array y for gemv" )

    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooJA.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( cooIA.size(), cooValues.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = cooIA.size();

    static LAMAKernel<COOKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = prefLoc;
    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<IndexType> rIA( cooIA, loc );
    ReadAccess<IndexType> rJA( cooJA, loc );
    ReadAccess<ValueType> rValues( cooValues, loc );

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<ValueType> rY( y, loc );
    WriteOnlyAccess<ValueType> wResult( result, loc, op == common::MatrixOp::NORMAL ? numRows : numColumns );

    normalGEMV[loc]( wResult.get(),
                     alpha, rX.get(), beta, rY.get(),
                     numRows, numColumns, numValues, rIA.get(), rJA.get(), rValues.get(), op );

    if ( async )
    {
        syncToken->pushRoutine( wResult.releaseDelayed() );
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rY.releaseDelayed() );
        syncToken->pushRoutine( rIA.releaseDelayed() );
        syncToken->pushRoutine( rJA.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

#define COOUTILS_SPECIFIER( ValueType )            \
                                                   \
    template void COOUtils::normalize(             \
            HArray<IndexType>&,                    \
            HArray<IndexType>&,                    \
            HArray<ValueType>&,                    \
            common::BinaryOp,                      \
            ContextPtr );                          \
                                                   \
    template void COOUtils::unique(                \
            HArray<IndexType>&,                    \
            HArray<IndexType>&,                    \
            HArray<ValueType>&,                    \
            common::BinaryOp,                      \
            ContextPtr );                          \
                                                   \
    template void COOUtils::sort(                  \
            HArray<IndexType>&,                    \
            HArray<IndexType>&,                    \
            HArray<ValueType>&,                    \
            ContextPtr );                          \
                                                   \
    template void COOUtils::getDiagonal(           \
            HArray<ValueType>&,                    \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            ContextPtr );                          \
                                                   \
    template void COOUtils::setDiagonal(           \
            HArray<ValueType>&,                    \
            const ValueType,                       \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            ContextPtr );                          \
                                                   \
    template void COOUtils::setDiagonalV(          \
            HArray<ValueType>&,                    \
            const HArray<ValueType>&,              \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            ContextPtr );                          \
                                                   \
    template void COOUtils::setRows(               \
            HArray<ValueType>&,                    \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            const common::BinaryOp op,             \
            ContextPtr );                          \
                                                   \
    template void COOUtils::setColumns(            \
            HArray<ValueType>&,                    \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            const common::BinaryOp op,             \
            ContextPtr );                          \
                                                   \
    template void COOUtils::binaryOp(              \
            HArray<IndexType>&,                    \
            HArray<IndexType>&,                    \
            HArray<ValueType>&,                    \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            const IndexType,                       \
            const IndexType,                       \
            const common::BinaryOp op,             \
            ContextPtr );                          \
                                                   \
    template void COOUtils::jacobi(                \
            HArray<ValueType>&,                    \
            const ValueType,                       \
            const HArray<ValueType>&,              \
            const HArray<ValueType>&,              \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            ContextPtr );                          \
                                                   \
    template void COOUtils::jacobiHalo(            \
            HArray<ValueType>&,                    \
            const ValueType,                       \
            const HArray<ValueType>&,              \
            const HArray<ValueType>&,              \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            ContextPtr );                          \
                                                   \
    template SyncToken* COOUtils::gemv(            \
            HArray<ValueType>&,                    \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            const common::MatrixOp,                \
            const ValueType,                       \
            const HArray<ValueType>&,              \
            const ValueType,                       \
            const HArray<ValueType>&,              \
            const bool,                            \
            ContextPtr );                          \

// selectComplexPart uses Math::real and Math::imag that is not defined for IndexType

SCAI_COMMON_LOOP( COOUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef COOUTILS_SPECIFIER

} /* end namespace sparsekernel */

} /* end namespace scai */
