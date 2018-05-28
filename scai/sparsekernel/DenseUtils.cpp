/**
 * @file DenseUtils.cpp
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
 * @brief Implementation and instantiation of utility methods for dense storage.
 * @author Thomas Brandes
 * @date 29.05.2018
 */

#include <scai/sparsekernel/DenseUtils.hpp>

#include <scai/sparsekernel/DenseKernelTrait.hpp>
#include <scai/sparsekernel/CSRUtils.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>
#include <scai/hmemo/HostReadAccess.hpp>

#include <scai/tracing.hpp>
#include <scai/common/macros/loop.hpp>

namespace scai
{

using namespace hmemo;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( DenseUtils::logger, "DenseUtils" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::convertDense2CSR(
    HArray<IndexType>& csrIA,
    HArray<IndexType>& csrJA,
    HArray<ValueType>& csrValues,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<ValueType>& denseValues,
    ContextPtr prefLoc )
{
    getSparseRowSizes( csrIA, numRows, numColumns, denseValues, prefLoc );

    SCAI_REGION( "Sparse.Dense.buildCSR" )

    IndexType numValues = HArrayUtils::scan1( csrIA, prefLoc );

    static LAMAKernel<DenseKernelTrait::getCSRValues<ValueType> > getCSRValues;

    ContextPtr loc = prefLoc;
    getCSRValues.getSupportedContext( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    WriteOnlyAccess<IndexType> wJA( csrJA, loc, numValues );
    WriteOnlyAccess<ValueType> wValues( csrValues, loc, numValues );

    ReadAccess<ValueType> rDense( denseValues, loc );

    SCAI_CONTEXT_ACCESS( loc )

    RealType<ValueType> eps = 0;

    getCSRValues[loc]( wJA.get(), wValues.get(), rIA.get(), numRows, numColumns,
                       rDense.get(), eps );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::getSparseRowSizes(
    HArray<IndexType>& csrSizes,
    const IndexType numRows,
    const IndexType numColumns,
    const hmemo::HArray<ValueType>& denseValues,
    ContextPtr prefLoc )
{
    SCAI_ASSERT_EQ_DEBUG( denseValues.size(), numRows * numColumns, 
                          "denseValues: size not " << numRows << " x " << numColumns )

    static LAMAKernel<DenseKernelTrait::getCSRSizes<ValueType> > getCSRSizes;

    ContextPtr loc = prefLoc;
    getCSRSizes.getSupportedContext( loc );
    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rDenseValues( denseValues, loc );

    // allocate sizes with one more element, so it can be used for offsets
    WriteOnlyAccess<IndexType> wSizes( csrSizes, loc, numRows + 1 );

    RealType<ValueType> eps = 0;

    getCSRSizes[loc]( wSizes.get(), numRows, numColumns, rDenseValues.get(), eps );
    wSizes.resize( numRows );
}

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void DenseUtils::convertCSR2Dense(
    HArray<ValueType>& denseValues,
    const IndexType numRows,
    const IndexType numColumns,
    const HArray<IndexType>& csrIA,
    const HArray<IndexType>& csrJA,
    const HArray<ValueType>& csrValues,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.Dense.setCSR" )

    static LAMAKernel<DenseKernelTrait::setCSRValues<ValueType> > setCSRValues;

    ContextPtr loc = prefLoc;
    setCSRValues.getSupportedContext( loc );

    SCAI_CONTEXT_ACCESS( loc );

    ReadAccess<IndexType> rIA( csrIA, loc );
    ReadAccess<IndexType> rJA( csrJA, loc );
    ReadAccess<ValueType> rValues( csrValues, loc );

    WriteOnlyAccess<ValueType> wDense( denseValues, loc, numRows * numColumns );

    setCSRValues[loc]( wDense.get(), numRows, numColumns, rIA.get(), rJA.get(), rValues.get() );
}

/* -------------------------------------------------------------------------- */

#define DENSE_UTILS_SPECIFIER( ValueType )           \
                                                     \
    template void DenseUtils::getSparseRowSizes(     \
        HArray<IndexType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::convertDense2CSR(      \
        HArray<IndexType>&,                          \
        HArray<IndexType>&,                          \
        HArray<ValueType>&,                          \
        const IndexType,                             \
        const IndexType,                             \
        const HArray<ValueType>&,                    \
        ContextPtr );                                \
                                                     \
    template void DenseUtils::convertCSR2Dense(      \
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
