/**
 * @file StencilUtils.hpp
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
 * @brief Utility functions for the STENCIL storage format.
 * @author Thomas Brandes
 * @date 01.10.2018
 */
#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/common/Grid.hpp>
#include <scai/common/Stencil.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/hmemo/Context.hpp>

namespace scai
{

namespace sparsekernel
{

/**
 *  This class provides a lot of utility functions regarding the stencil storage format.
 *
 *  A stencil storage is described by:
 *
 *   - the grid to which it is applied to
 *   - the stencil that specifies the linear operation
 *
 *  When generating a stencil storage some heterogeneous arrays are set up that contain
 *  the relevant infomation.
 *
 *   - gridArray contains grid width, grid distances
 *   - stencilGeometry contains stencil positions and distances
 *   - stencilValues contains the values
 */
class COMMON_DLL_IMPORTEXPORT StencilUtils
{
public:

    /** 
     *  @brief Method to setup heterogeneous arrays that contain all relevant grid and stencil info.
     *
     *  @param[in] grid (nDims) contains grid sizes for the actual grid to which stencil is applied
     *  @param[in] stencil contains the stencil positions and values (nPoints)
     *  @param[out] gridInfo contains gridSizes (nDims), gridDistances(nDims), gridBorders(2*nDims), gridStencilWidth(2*nDims)
     *  @param[out] stencilInfo contains stencilPositions (nDims*nPoints), stencilOffsets(nPoints)
     *  @param[out] stencilValues contains stencil values[nPoints]
     *
     *  The heterogeneous arrays are required to call stencil kernels on arbitrary devices and to avoid recomputation of
     *  data that is reused, e.g. for multiple matrix-vector multiplications.
     */
    template<typename ValueType>
    static void setup(
        hmemo::HArray<IndexType>& gridInfo,
        hmemo::HArray<int>& stencilInfo,
        hmemo::HArray<ValueType>& stencilValues,
        const common::Grid& grid,
        const common::Stencil<ValueType>& stencil );

    template<typename ValueType>
    static tasking::SyncToken* gemv(
        hmemo::HArray<ValueType>& result,
        const ValueType alpha,
        const hmemo::HArray<ValueType>& x,
        const ValueType beta,
        const hmemo::HArray<ValueType>& y,
        const IndexType gridSize,
        const IndexType nDim,
        const IndexType gridSizes[],
        const IndexType nPoints,
        const hmemo::HArray<IndexType>& gridInfo,
        const hmemo::HArray<int>& stencilInfo,
        const hmemo::HArray<ValueType>& stencilValues,
        bool async,
        hmemo::ContextPtr prefLoc );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

/* -------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */
