/**
 * @file StencilUtils.cpp
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
 * @brief Implementation and instantiation of utility methods for Stencil storage.
 * @author Thomas Brandes
 * @date 01.10.2018
 */

#include <scai/sparsekernel/StencilUtils.hpp>

#include <scai/sparsekernel/StencilKernelTrait.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>
#include <scai/utilskernel/LAMAKernel.hpp>

#include <scai/common/macros/loop.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using namespace hmemo;

using utilskernel::LAMAKernel;
using utilskernel::HArrayUtils;
using tasking::SyncToken;
using common::Grid;
using common::Stencil;

namespace sparsekernel
{

SCAI_LOG_DEF_LOGGER( StencilUtils::logger, "StencilUtils" )

/* -------------------------------------------------------------------------- */

template<typename ValueType>
void StencilUtils::setup(
    HArray<IndexType>& gridInfo,
    HArray<int>& stencilInfo,
    HArray<ValueType>& stencilValues,
    const Grid& grid,
    const Stencil<ValueType>& stencil )
{
    SCAI_ASSERT_EQ_ERROR( grid.nDims(), stencil.nDims(), "Serious dimension mismatch grid/stencil" )

    IndexType nDims = grid.nDims();

    // gridInfo[6*nDims] = gridSizes[nDims] | gridDistance[nDims] | gridBorders[2*nDims] | gridStencilWidth[2*nDims]

    // Note: setup is done completely on the host, arrays might be copied to other context later

    auto wGridInfo = hostWriteOnlyAccess( gridInfo, 6 * nDims );
   
    const IndexType* pGridSizes = grid.sizes();
    const Grid::BorderType* pGridBorders = grid.borders();

    for ( IndexType i = 0; i < nDims; ++i )
    {
        wGridInfo[i] = pGridSizes[i];
        wGridInfo[ 2 * nDims + 2 * i] = static_cast<IndexType>( pGridBorders[2 * i] );
        wGridInfo[ 2 * nDims + 2 * i + 1] = static_cast<IndexType>( pGridBorders[2 * i + 1] );
    }

    grid.getDistances( wGridInfo.get() + nDims );
    stencil.getWidth( wGridInfo.get() + 4 * nDims );

    const IndexType nPoints = stencil.nPoints();

    // stencilInfo[nDims*nPoints + nPoints] = stenilData[nDims*nPoints] | stencilOffsets[nPoints]

    auto wStencilInfoPure = hostWriteOnlyAccess( stencilInfo, nPoints * nDims + nPoints );

    // current workaround required as there might be no HArray<int>

    int* wStencilInfo = reinterpret_cast<int*>( wStencilInfoPure.get() );

    const int* pStencilPositions = stencil.positions();
   
    for ( IndexType i = 0; i < nPoints * nDims; ++ i )
    {
        wStencilInfo[i] = pStencilPositions[i];
    }

    stencil.getLinearOffsets( wStencilInfo + nPoints * nDims, wGridInfo.get() + nDims );

    auto wStencilValues = hostWriteOnlyAccess( stencilValues, nPoints );

    const ValueType* pStencilValues = stencil.values();

    for ( IndexType i = 0; i < nPoints; ++ i )
    {
        wStencilValues[i] = pStencilValues[i];
    }
}

template<typename ValueType>
SyncToken* StencilUtils::gemv(
    HArray<ValueType>& result,
    const ValueType alpha,
    const HArray<ValueType>& x,
    const ValueType beta,
    const HArray<ValueType>& y,
    const IndexType gridSize,
    const IndexType nDims,
    const IndexType gridSizes[],
    const IndexType nPoints,
    const HArray<IndexType>& gridInfo,
    const HArray<int>& stencilInfo,
    const HArray<ValueType>& stencilValues,
    bool async,
    ContextPtr prefLoc )
{
    SCAI_REGION( "Sparse.Stencil.gemv" )

    SCAI_ASSERT_EQ_ERROR( x.size(), gridSize, "serious size mismatch" )
    SCAI_ASSERT_EQ_ERROR( gridInfo.size(), 6 * nDims, "serious mismatch"  );
    SCAI_ASSERT_EQ_ERROR( stencilInfo.size(), nPoints * ( nDims + 1 ), 
                             "serious mismatch for " << nDims << "D" << nPoints << "P stencil" );
    SCAI_ASSERT_EQ_ERROR( stencilValues.size(), nPoints, "serious mismatch" );

    static LAMAKernel<StencilKernelTrait::normalGEMV<ValueType> > normalGEMV;

    ContextPtr loc = prefLoc;

    normalGEMV.getSupportedContext( loc );

    std::unique_ptr<SyncToken> syncToken;

    if ( async )
    {
        syncToken.reset( loc->getSyncToken() );
    }

    SCAI_ASYNCHRONOUS( syncToken.get() );

    SCAI_CONTEXT_ACCESS( loc )

    ReadAccess<ValueType> rX( x, loc );
    ReadAccess<IndexType> rGridInfo( gridInfo, loc );
    ReadAccess<int> rStencilInfo( stencilInfo, loc );
    ReadAccess<ValueType> rStencilValues( stencilValues, loc );

    // get the different grid info data as context pointers

    const IndexType* dGridSizes = rGridInfo.get();
    const IndexType* dDistances = dGridSizes + nDims;
    const IndexType* dBorders = dDistances + nDims;
    const IndexType* dStencilWidth = dBorders + 2 * nDims;

    // get the different stencil info data as context pointers

    const int* dStencilPositions = rStencilInfo.get();
    const int* dStencilOffsets   = dStencilPositions + nDims * nPoints;

    if ( beta != ValueType( 0 ) )
    {
        SCAI_ASSERT_EQ_ERROR( y.size(), gridSize, "y has illegal size" )

        ReadAccess<ValueType> rY( y, loc );
        WriteOnlyAccess<ValueType> wResult( result, loc, gridSize );  // result might be aliased to y

        SCAI_LOG_INFO( logger, "call kernel normalGEMV( beta = " << beta << ", y[" << gridSize << "] = " << rY.get() << " ) on " << *loc )

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), beta, rY.get(), 
                         nDims, gridSizes, dGridSizes, dDistances, dBorders, dStencilWidth,
                         nPoints, dStencilPositions, rStencilValues.get(), dStencilOffsets );

        if ( async )
        {
            syncToken->pushRoutine( rY.releaseDelayed() );
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }
    else
    {
        // beta == 0, do not access y at all

        WriteOnlyAccess<ValueType> wResult( result, loc, gridSize );

        SCAI_LOG_INFO( logger, "call kernel normalGEMV( beta is 0 ) on " << *loc )

        normalGEMV[loc]( wResult.get(), alpha, rX.get(), ValueType( 0 ), NULL,
                         nDims, gridSizes, dGridSizes, dDistances, dBorders, dStencilWidth,
                         nPoints, dStencilPositions, rStencilValues.get(), dStencilOffsets );

        if ( async )
        {
            syncToken->pushRoutine( wResult.releaseDelayed() );
        }
    }

    if ( async )
    {
        syncToken->pushRoutine( rX.releaseDelayed() );
        syncToken->pushRoutine( rGridInfo.releaseDelayed() );
        syncToken->pushRoutine( rStencilInfo.releaseDelayed() );
        syncToken->pushRoutine( rStencilValues.releaseDelayed() );
    }

    return syncToken.release();
}

/* -------------------------------------------------------------------------- */

#define STENCIL_UTILS_SPECIFIER( ValueType )         \
                                                     \
    template void StencilUtils::setup(               \
        HArray<IndexType>&,                          \
        HArray<int>&,                                \
        HArray<ValueType>&,                          \
        const Grid&,                                 \
        const Stencil<ValueType>& );                 \
                                                     \
    template SyncToken* StencilUtils::gemv(          \
        HArray<ValueType>&,                          \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const ValueType,                             \
        const HArray<ValueType>&,                    \
        const IndexType,                             \
        const IndexType,                             \
        const IndexType[],                           \
        const IndexType,                             \
        const HArray<IndexType>&,                    \
        const HArray<int>&,                          \
        const HArray<ValueType>&,                    \
        const bool,                                  \
        ContextPtr );                                \
                                                     \

SCAI_COMMON_LOOP( STENCIL_UTILS_SPECIFIER, SCAI_NUMERIC_TYPES_HOST )

#undef STENCIL_UTILS_SPECIFIER

} /* end namespace utilskernel */

} /* end namespace scai */
