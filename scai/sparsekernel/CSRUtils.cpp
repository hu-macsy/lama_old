/**
 * @file CSRUtils.cpp
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
 * @brief Implementation and instantion of CSR utility methods.
 * @author Thomas Brandes
 * @date 14.02.2018
 */

#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/sparsekernel/CSRKernelTrait.hpp>

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

SCAI_LOG_DEF_LOGGER( CSRUtils::logger, "CSRUtils" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void CSRUtils::sort(
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values,
    const hmemo::HArray<IndexType>& ia,
    const IndexType numColumns,
    const bool diagonalFlag,
    const ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_ERROR( values.size(), ja.size(), "serious size mismatch" )

    const IndexType numValues = values.size();   // #non-zero entries
    const IndexType numRows = ia.size() - 1;

    static LAMAKernel<CSRKernelTrait::sortRowElements<ValueType> > sortRowElements;

    ContextPtr loc = prefLoc;
    sortRowElements.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    WriteAccess<IndexType> wJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    sortRowElements[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numColumns, numValues, diagonalFlag );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
IndexType CSRUtils::setDiagonalFirst(
    hmemo::HArray<IndexType>& ja,
    hmemo::HArray<ValueType>& values,
    const hmemo::HArray<IndexType>& ia,
    const IndexType numColumns,
    hmemo::ContextPtr prefLoc )
{
    const IndexType numRows = ia.size() - 1;

    static LAMAKernel<CSRKernelTrait::setDiagonalFirst<ValueType> > setDiagonal;

    ContextPtr loc = prefLoc;
    setDiagonal.getSupportedContext( loc );

    SCAI_LOG_INFO( logger, "setDiagonalFirst on CSR data ( " << numRows << " x " << numColumns
                     << " ), called on " << *loc << ", preferred was " << *prefLoc )

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<IndexType> rIA( ia, loc );
    WriteAccess<IndexType> wJA( ja, loc );
    WriteAccess<ValueType> wValues( values, loc );

    IndexType numDiagonals = std::min( numRows, numColumns );

    return setDiagonal[loc]( wJA.get(), wValues.get(), numDiagonals, rIA.get() );
}

/* -------------------------------------------------------------------------- */

#define CSRUTILS_SPECIFIER( ValueType )              \
    template void CSRUtils::sort(                    \
            hmemo::HArray<IndexType>&,               \
            hmemo::HArray<ValueType>&,               \
            const hmemo::HArray<IndexType>&,         \
            const IndexType,                         \
            const bool,                              \
            ContextPtr );                            \
    template IndexType CSRUtils::setDiagonalFirst(   \
            hmemo::HArray<IndexType>&,               \
            hmemo::HArray<ValueType>&,               \
            const hmemo::HArray<IndexType>&,         \
            const IndexType,                         \
            ContextPtr );                            \

SCAI_COMMON_LOOP( CSRUTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef CSRUTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
