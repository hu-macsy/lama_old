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

#include <scai/lama/storage/COOUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <scai/common/macros/loop.hpp>
#include <algorithm>

namespace scai
{

using namespace hmemo;

namespace lama
{

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void COOUtils::sort(
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values )
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
        cmp( const hmemo::HArray<IndexType>& ia, const hmemo::HArray<IndexType>& ja )
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
    hmemo::HArray<IndexType>& ia,
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values,
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

#define COOUTILS_SPECIFIER( ValueType )            \
    template void COOUtils::unique(                \
            hmemo::HArray<IndexType>&,             \
            hmemo::HArray<IndexType>&,             \
            hmemo::HArray<ValueType>&,             \
            common::BinaryOp );                    \
    template void COOUtils::sort(                  \
            hmemo::HArray<IndexType>&,             \
            hmemo::HArray<IndexType>&,             \
            hmemo::HArray<ValueType>& );           \

// selectComplexPart uses Math::real and Math::imag that is not defined for IndexType

SCAI_COMMON_LOOP( COOUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef COOUTILS_SPECIFIER

} /* end namespace lama */

} /* end namespace scai */
