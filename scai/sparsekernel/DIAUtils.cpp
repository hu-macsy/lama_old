/**
 * @file DIAUtils.cpp
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
 * @brief Implementation and instantiation of utility methods for DIA storage.
 * @author Thomas Brandes
 * @date 29.05.2018
 */

#include <scai/sparsekernel/DIAUtils.hpp>

#include <scai/sparsekernel/DIAKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>

#include <memory>

namespace scai
{

using namespace hmemo;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using tasking::SyncToken;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( DIAUtils::logger, "DIAUtils" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIAUtils::convertDIA2CSR(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& diaOffset,
    const HArray<ValueType>& diaValues,
    ContextPtr prefLoc )
{
    getCSRSizes( csrIA, numRows, numColumns, diaOffset, diaValues, prefLoc );

    const IndexType numDiagonals = diaOffset.size();

    IndexType numValues = HArrayUtils::scan1( csrIA, prefLoc );

    static LAMAKernel<DIAKernelTrait::getCSRValues<ValueType> > getCSRValues;
    ContextPtr loc = prefLoc;
    getCSRValues.getSupportedContext( loc );

    ReadAccess<ValueType> rValues( diaValues, loc );
    ReadAccess<IndexType> rOffset( diaOffset, loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    WriteOnlyAccess<IndexType> wJA( csrJA, loc, numValues );
    WriteOnlyAccess<ValueType> wValues( csrValues, loc, numValues );

    SCAI_CONTEXT_ACCESS( loc )

    getCSRValues[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numColumns, numDiagonals,
                       rOffset.get(), rValues.get() );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIAUtils::getCSRSizes(
    HArray<IndexType>& csrSizes,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<OffsetType>& diaOffset,
    const hmemo::HArray<ValueType>& diaValues,
    ContextPtr prefLoc )
{
    static LAMAKernel<DIAKernelTrait::getCSRSizes<ValueType> > getCSRSizes;

    ContextPtr loc = prefLoc;
    getCSRSizes.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    IndexType numDiagonals = diaOffset.size();

    ReadAccess<ValueType> rValues( diaValues, loc );
    ReadAccess<IndexType> rOffset( diaOffset, loc );

    WriteOnlyAccess<IndexType> wIA( csrSizes, loc, numRows + 1 );

    getCSRSizes[loc]( wIA.get(), numRows, numColumns, numDiagonals, rOffset.get(), rValues.get() );

    wIA.resize( numRows );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIAUtils::getDIAOffset(
    HArray<IndexType>& diaOffset,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    ContextPtr )
{
    // build a set of all used lower and upper diagonals

    IndexType maxNumDiagonals = common::Math::max( numRows, numColumns );
    IndexType allNumDiagonals = 2 * maxNumDiagonals - 1;
    IndexType mainDiagonal = maxNumDiagonals - 1;

    std::unique_ptr<bool[]> diagonalUsed( new bool[allNumDiagonals] );

    for ( IndexType i = 0; i < allNumDiagonals; i++ )
    {
        diagonalUsed[i] = false;
    }

    auto rIA = hostReadAccess( csrIA );
    auto rJA = hostReadAccess( csrJA );
    auto rValues = hostReadAccess( csrValues );

    #pragma omp parallel for 

    for ( IndexType i = 0; i < numRows; ++i )
    {
        for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; jj++ )
        {
            IndexType j = rJA[jj]; // column used

            bool& flag = diagonalUsed[ mainDiagonal + ( j - i ) ];

            if ( !flag )
            {
                flag = true;    
            }
        }
    }

    IndexType numDiagonals = 0;

    for ( IndexType i = 0; i < allNumDiagonals; i++ )
    {
        if ( diagonalUsed[i] )
        {   
            numDiagonals++;
        }
    }

    if ( numDiagonals == 0 )
    {
        diaOffset.clear();
        return;
    }

    WriteOnlyAccess<IndexType> wOffset( diaOffset, numDiagonals );

    numDiagonals = 0;

    for ( OffsetType i = 0; i < allNumDiagonals; ++i )
    {
        if ( diagonalUsed[i] )
        {
            wOffset[numDiagonals++] = i - mainDiagonal;
        }
    }
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIAUtils::convertCSR2DIA(
    HArray<IndexType>& diaOffset,
    HArray<ValueType>& diaValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    ContextPtr prefLoc )
{
    SCAI_LOG_INFO( logger, "convert CSR " << numRows << " x " << numColumns
                     << ", ia = " << csrIA << ", ja = " << csrJA << ", values = " << csrValues )

    if ( numRows == 0 || numColumns == 0 || csrJA.size() == 0 )
    {
        // no entries, no diagonals, avoid bad allocations

        diaOffset.clear();
        diaValues.clear();
        return;
    }

    getDIAOffset( diaOffset, numRows, numColumns, csrIA, csrJA, csrValues, prefLoc );

    const IndexType numDiagonals = diaOffset.size();

    ValueType ZERO = 0;

    auto rIA = hostReadAccess( csrIA );
    auto rJA = hostReadAccess( csrJA );
    auto rValues = hostReadAccess( csrValues );

    auto wValues = hostWriteOnlyAccess( diaValues, numRows * numDiagonals );
    auto rOffset = hostReadAccess( diaOffset );

    #pragma omp parallel for 
    for ( IndexType i = 0; i < numRows; i++ )
    {
        for ( IndexType d = 0; d < numDiagonals; d++ )
        {
            ValueType& addrValue = wValues[ d * numRows + i ];

            addrValue = ZERO;

            IndexType j = i + rOffset[d];

            if ( !common::Utils::validIndex( j, numColumns ) )
            {
                continue;
            }

            for ( IndexType jj = rIA[i]; jj < rIA[i + 1]; ++jj )
            {
                if ( rJA[jj] == j )
                {
                    addrValue = static_cast<ValueType>( rValues[jj] );
                    break;
                }
            }
        }
    }
}

/* -------------------------------------------------------------------------- */

IndexType DIAUtils::getValuePos(
    const IndexType i,
    const IndexType j,
    const IndexType numRows,
    const IndexType,
    const hmemo::HArray<OffsetType>& diaOffset,
    hmemo::ContextPtr prefLoc )
{
    static LAMAKernel<DIAKernelTrait::getValuePos> getValuePos;

    ContextPtr loc = prefLoc;
    getValuePos.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rOffset( diaOffset, loc );

    IndexType numDiagonals = diaOffset.size();
    IndexType pos = getValuePos[loc]( i, j, numRows, rOffset.get(), numDiagonals );

    return pos;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIAUtils::jacobi(
    hmemo::HArray<ValueType>& solution,
    const ValueType omega,
    const hmemo::HArray<ValueType>& oldSolution,
    const hmemo::HArray<ValueType>& rhs,
    const IndexType n,
    const hmemo::HArray<OffsetType>& diaOffset,
    const hmemo::HArray<ValueType>& diaValues,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( n, oldSolution.size(), "size mismatch" )
    SCAI_ASSERT_EQ_DEBUG( n, rhs.size(), "size mismatch" )

    static LAMAKernel<DIAKernelTrait::jacobi<ValueType> > jacobi;

    const IndexType numDiagonals = diaOffset.size();

    ContextPtr loc = prefLoc;
    jacobi.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> rOffset( diaOffset, loc );
    ReadAccess<ValueType> rValues( diaValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rRhs( rhs, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, n );
    jacobi[loc]( wSolution.get(), n, numDiagonals, rOffset.get(), rValues.get(),
                 rOldSolution.get(), rRhs.get(), omega );
}

/* -------------------------------------------------------------------------- */

void DIAUtils::getRowPositions(
    hmemo::HArray<IndexType>& indexes,
    hmemo::HArray<IndexType>& positions,
    const IndexType i,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<OffsetType>& diaOffset,
    hmemo::ContextPtr )
{
    IndexType numDiagonals = diaOffset.size();

    // numDiagonals is also maximal number of non-zero entries in a row 

    WriteOnlyAccess<IndexType> wIndexes( indexes, numDiagonals );
    WriteOnlyAccess<IndexType> wPositions( positions, numDiagonals );

    const ReadAccess<OffsetType> rOffset( diaOffset );

    IndexType cnt = 0;  // count the real number of entries, might be lt numDiagonals 

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        IndexType j = i + rOffset[d];

        if ( common::Utils::validIndex( j, numColumns ) )
        {
            wIndexes[cnt] = j;
            wPositions[cnt] = d * numRows + i;
            ++cnt;
        }
    }

    wIndexes.resize( cnt );
    wPositions.resize( cnt );
}

/* -------------------------------------------------------------------------- */

void DIAUtils::getColPositions(
    hmemo::HArray<IndexType>& indexes,
    hmemo::HArray<IndexType>& positions,
    const IndexType j,
    const IndexType numRows,
    const IndexType ,
    const hmemo::HArray<OffsetType>& diaOffset,
    hmemo::ContextPtr )
{
    IndexType numDiagonals = diaOffset.size();
    IndexType sizeDiagonal = numRows;

    // numDiagonals is also maximal number of non-zero entries in a column

    WriteOnlyAccess<IndexType> wIndexes( indexes, numDiagonals );
    WriteOnlyAccess<IndexType> wPositions( positions, numDiagonals );
    
    const ReadAccess<IndexType> rOffset( diaOffset );

    IndexType cnt = 0;  // count the real number of entries, might be lt numDiagonals 

    for ( IndexType d = 0; d < numDiagonals; ++d )
    {
        IndexType i = j - rOffset[d];

        if ( common::Utils::validIndex( i, numRows ) )
        {
            wIndexes[cnt] = i;
            wPositions[cnt] = d * sizeDiagonal + i;
            ++cnt;
        }
    }

    wIndexes.resize( cnt );
    wPositions.resize( cnt );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DIAUtils::jacobiHalo(
    hmemo::HArray<ValueType>& solution,
    const ValueType omega,
    const hmemo::HArray<ValueType>& diagonal,
    const hmemo::HArray<ValueType>& oldSolution,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<OffsetType>& diaOffset,
    const hmemo::HArray<ValueType>& diaValues,
    hmemo::ContextPtr prefLoc )
{
    static LAMAKernel<DIAKernelTrait::jacobiHalo<ValueType> > jacobiHalo;

    const IndexType numDiagonals = diaOffset.size();

    ContextPtr loc = prefLoc;
    jacobiHalo.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )
    ReadAccess<IndexType> rOffset( diaOffset, loc );
    ReadAccess<ValueType> rValues( diaValues, loc );
    ReadAccess<ValueType> rOldSolution( oldSolution, loc );
    ReadAccess<ValueType> rDiagonal( diagonal, loc );
    WriteOnlyAccess<ValueType> wSolution( solution, loc, numRows );

    jacobiHalo[loc]( wSolution.get(), rDiagonal.get(),
                     numRows, numColumns, numDiagonals, 
                     rOffset.get(), rValues.get(),
                     rOldSolution.get(), omega );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
RealType<ValueType> DIAUtils::maxNorm(
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<OffsetType>& diaOffset,
    const hmemo::HArray<ValueType>& diaValues,
    hmemo::ContextPtr prefLoc )
{
    IndexType numDiagonals = diaOffset.size();

    static LAMAKernel<DIAKernelTrait::absMaxVal<ValueType> > absMaxVal;

    ContextPtr loc = prefLoc;

    absMaxVal.getSupportedContext( loc );

    ReadAccess<IndexType> rOffsets( diaOffset, loc );
    ReadAccess<ValueType> rValues( diaValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    RealType<ValueType> maxval = 
        absMaxVal[loc]( numRows, numColumns, numDiagonals, rOffsets.get(), rValues.get() );

    return maxval;
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
SyncToken* DIAUtils::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<OffsetType>& diaOffset,
    const HArray<ValueType>& diaValues,
    const common::MatrixOp op,
    bool async,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Storage.DIA.gemv" )

    const IndexType nSource = common::isTranspose( op ) ? numRows : numColumns;
    const IndexType nTarget = common::isTranspose( op ) ? numColumns : numRows;

    IndexType numDiagonals = diaOffset.size();

    SCAI_ASSERT_EQUAL_ERROR( x.size(), nSource )

    static LAMAKernel<DIAKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = prefLoc;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rOffset( diaOffset, loc );
    ReadAccess<ValueType> rValues( diaValues, loc );
    ReadAccess<ValueType> rX( x, loc );

    if ( beta != ValueType( 0 ) )
    {
        SCAI_ASSERT_EQ_ERROR( y.size(), nTarget, "y has illegal size" )

        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );  // result might be aliased to y

        SCAI_LOG_INFO( logger, "call kernel normalGEMV( beta = " << beta << ", y[" << nTarget << "] = " << rY.get() << " ) on " << *loc )

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(),
                         numRows, numColumns, numDiagonals,
                         rOffset.get(), rValues.get(), op );

        if ( async )
        {
            syncToken->pushRoutine( rY.releaseDelayed() );
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }
    else
    {
        // do not access y at all

        WriteOnlyAccess<ValueType> wResult( result, loc, nTarget );

        SCAI_LOG_INFO( logger, "call kernel normalGEMV( beta is 0 ) on " << *loc )

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL,
                         numRows, numColumns, numDiagonals,
                         rOffset.get(), rValues.get(), op );

        if ( async )
        {
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }

    if ( async )
    {
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rValues.releaseDelayed() );
        syncToken->pushRoutine( rOffset.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

#define DENSE_UTILS_SPECIFIER( ValueType )           \
                                                     \
    template void DIAUtils::convertDIA2CSR(          \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DIAUtils::convertCSR2DIA(          \
        HArray<IndexType>&,                          \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<IndexType>&,                    \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DIAUtils::jacobi(                  \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const HArray<OffsetType>&,                   \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DIAUtils::jacobiHalo(              \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<OffsetType>&,                   \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template RealType<ValueType> DIAUtils::maxNorm(  \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<OffsetType>&,                   \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template SyncToken* DIAUtils::gemv(              \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<OffsetType>&,                   \
        const HArray<ValueType>&,                    \
        const common::MatrixOp op,                   \
        const bool,                                  \
        ContextPtr );                                \

SCAI_COMMON_LOOP( DENSE_UTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef DENSE_UTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
