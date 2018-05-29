/**
 * @file DIAUtils.cpp
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
        for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; jj++ )
        {
            IndexType j = csrJA[jj]; // column used

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
    SCAI_LOG_ERROR( logger, "convert CSR " << numRows << " x " << numColumns
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

            for ( IndexType jj = csrIA[i]; jj < csrIA[i + 1]; ++jj )
            {
                if ( csrJA[jj] == j )
                {
                    addrValue = static_cast<ValueType>( csrValues[jj] );
                    break;
                }
            }
        }
    }
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

SCAI_COMMON_LOOP( DENSE_UTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef DENSE_UTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
