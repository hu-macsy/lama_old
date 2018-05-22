/**
 * @file COOUtils.cpp
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
 * @brief Utility functions for COO data
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/COOUtils.hpp>
#include <scai/sparsekernel/COOKernelTrait.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <scai/common/macros/loop.hpp>
#include <algorithm>

namespace scai
{

using namespace hmemo;

using utilskernel::LAMAKernel;

namespace sparsekernel
{

/* -------------------------------------------------------------------------- */

void COOUtils::convertCOO2CSR(
    hmemo::HArray<IndexType>& csrIA,
    const HArray<IndexType>& cooIA,
    const IndexType numRows,
    ContextPtr )
{
    const IndexType nnz = cooIA.size();

    ContextPtr loc = Context::getHostPtr();

    WriteOnlyAccess<IndexType> wOffsets( csrIA, loc, numRows + 1 );
    ReadAccess<IndexType> rIA( cooIA, loc );
   
    wOffsets[0] = 0;

    for ( IndexType i = 0; i < numRows; ++i )
    {
        IndexType offs = wOffsets[i];

        while ( offs < nnz && rIA[offs] == i )
        {
            offs++;
        }

        wOffsets[i + 1] = offs;
    }
}

/* -------------------------------------------------------------------------- */

void COOUtils::convertCSR2COO(
    hmemo::HArray<IndexType>& cooIA,
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

    auto rIA = hostReadAccess( ia );
    auto rJA = hostReadAccess( ja );

    bool sorted = true;

    #pragma parallel for

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
    HArray<ValueType>& values )
{
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
    common::BinaryOp op )

{
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
void COOUtils::getDiagonal(
    hmemo::HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    const hmemo::HArray<ValueType>& values,
    hmemo::ContextPtr prefLoc )
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
    hmemo::HArray<ValueType>& values,
    const hmemo::HArray<ValueType>& diagonal,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( ia.size(), values.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = ia.size();
    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<COOKernelTrait::setDiagonalV<ValueType>> kSetDiagonalV;

    // choose Context where all kernel routines are available

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
    hmemo::HArray<ValueType>& values,
    const ValueType diagonalValue,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    hmemo::ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( ia.size(), ja.size(), "serious size mismatch for COO arrays" )
    SCAI_ASSERT_EQ_ERROR( ia.size(), values.size(), "serious size mismatch for COO arrays" )

    const IndexType numValues = ia.size();
    const IndexType numDiagonals = common::Math::min( numRows, numColumns );

    static LAMAKernel<COOKernelTrait::setDiagonal<ValueType>> kSetDiagonal;

    // choose Context where all kernel routines are available

    ContextPtr loc = prefLoc;
    kSetDiagonal.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( ia, loc );
    ReadAccess<IndexType> rJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    SCAI_CONTEXT_ACCESS( loc )

    kSetDiagonal[loc]( wValues.get(), diagonalValue, numDiagonals, rIA.get(), rJA.get(), numValues );
}

/* -------------------------------------------------------------------------- */

bool COOUtils::hasDiagonalProperty(
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<IndexType>& ia,
    const hmemo::HArray<IndexType>& ja,
    hmemo::ContextPtr prefLoc )
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

#define COOUTILS_SPECIFIER( ValueType )            \
    template void COOUtils::unique(                \
            HArray<IndexType>&,                    \
            HArray<IndexType>&,                    \
            HArray<ValueType>&,                    \
            common::BinaryOp );                    \
    template void COOUtils::sort(                  \
            HArray<IndexType>&,                    \
            HArray<IndexType>&,                    \
            HArray<ValueType>& );                  \
    template void COOUtils::getDiagonal(           \
            HArray<ValueType>&,                    \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            const HArray<ValueType>&,              \
            ContextPtr );                          \
    template void COOUtils::setDiagonal(           \
            HArray<ValueType>&,                    \
            const ValueType,                       \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            ContextPtr );                          \
    template void COOUtils::setDiagonalV(          \
            HArray<ValueType>&,                    \
            const HArray<ValueType>&,              \
            const IndexType,                       \
            const IndexType,                       \
            const HArray<IndexType>&,              \
            const HArray<IndexType>&,              \
            ContextPtr );                          \

// selectComplexPart uses Math::real and Math::imag that is not defined for IndexType

SCAI_COMMON_LOOP( COOUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef COOUTILS_SPECIFIER

} /* end namespace sparsekernel */

} /* end namespace scai */
