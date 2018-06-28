/**
 * @file StencilKernelTrait.hpp
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
 * @brief Definition of kernel interfaces used for stencil operations.
 * @author Thomas Brandes
 * @date 27.04.2017
 */

#pragma once

#include <scai/common/Grid.hpp>

namespace scai
{
namespace sparsekernel
{

struct StencilKernelTrait
{
    struct stencilLocalSizes
    {
        /** Function that computes local csr sizes (num entries per point) for a stencil matrix.
         *
         *  @param[out] sizes is array with entry for each grid point how many stencil points are available
         *  @param[in] nDims is the number of dimensions for the grid
         *  @param[in] gridSizes contains the dimensions of the grid 
         *  @param[in] gridDistances contains the distance between two neighbored points for each dim
         *  @param[in] gridBorders contains for each dim two entries to specify kind of border
         *  @param[in] nPoints number of stencil points
         *  @param[in] stencilNodes contins nDims * nPoints direction values for the stencil points
         *
         *  The grid is assumed to be linearized for entries in sizes where the kind of linearization is 
         *  given by the gridDistances, e.g. [ numCols, 1 ] for row-major, [1, numRows ] for col-major.
         */
        typedef void ( *FuncType )(
            IndexType sizes[],
            const IndexType nDims,
            const IndexType gridSizes[],
            const IndexType gridDistances[],
            const common::Grid::BorderType gridBorders[],
            const IndexType nPoints,
            const int stencilNodes[] );

        static const char* getId()
        {
            return "Stencil.localSizes";
        }
    };

    template<typename ValueType>
    struct stencilLocalCSR
    {
        /** Function that builds local csr data for a stencil matrix.
         *
         *  @param[out] csrJA  contains the column indexes for the csr matrix
         *  @param[out] csrValues  contains the non-zero values for the csr matrix corresponding to csrJA
         *  @param[in] csrIA offsets array, csrIA[i] gives position where values start belonging to grid point i
         *  @param[in] nDims is the number of dimensions for the grid
         *  @param[in] gridSizes contains the dimensions of the grid 
         *  @param[in] gridDistances contains the distance between two neighbored points for each dim
         *  @param[in] gridBorders contains for each dim two entries to specify kind of border
         *  @param[in] nPoints number of stencil points
         *  @param[in] stencilNodes contins nDims * nPoints direction values for the stencil points
         *  @param[in] stencilVal contains the scale value for each stencil point
         *  @param[in] stencilOffset array with offset for each stencil point in linearized grid
         *
         *  csrIA is the offset array that is a scan of the sizes array computed by stencilSizes.
         */
        typedef void ( *FuncType )(
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType nDims,
            const IndexType gridSizes[],
            const IndexType gridDistances[],
            const common::Grid::BorderType gridBorders[],
            const IndexType nPoints,
            const int stencilNodes[],
            const ValueType stencilVal[],
            const int stencilOffset[] );
    
        static const char* getId()
        {
            return "Stencil.localCSR";
        }
    };

    template<typename ValueType>
    struct stencilHaloCSR
    {
        /** Function that build halo csr data for a stencil matrix.
         *
         *  @param[out] csrJA  contains the 'global' column indexes for the csr matrix
         *  @param[out] csrValues  contains the non-zero values for the csr matrix corresponding to csrJA
         *  @param[in] csrIA offsets array, csrIA[i] gives position where values start belonging to grid point i
         *  @param[in] nDims is the number of dimensions for the grid
         *  @param[in] localGridSizes contains the dimensions of the grid 
         *  @param[in] localGridDistances contains the distance between two neighbored points for each dim
         *  @param[in] localLB where the local grid is positioned in the global grid
         *  @param[in] globalGridSizes contains the dimensions of the global grid 
         *  @param[in] globalGridDistances contains the distance between two neighbored points for each dim
         *  @param[in] globalGridBorders contains for each dimension two entries specifying how to deal with border
         *  @param[in] nPoints number of stencil points
         *  @param[in] stencilNodes contins nDims * nPoints direction values for the stencil points
         *  @param[in] stencilVal contains the scale value for each stencil point
         *  @param[in] stencilOffset array with offset for each stencil point in linearized 'global' grid
         *
         *  csrIA is the offset array that is a scan of the sizes array computed by stencilSizes.
         */
        typedef void ( *FuncType )(
            IndexType csrJA[],
            ValueType csrValues[],
            const IndexType csrIA[],
            const IndexType nDims,
            const IndexType localGridSizes[],
            const IndexType localGridDistances[],
            const IndexType localLB[],
            const IndexType globalGridSizes[],
            const IndexType globalGridDistances[],
            const common::Grid::BorderType globalGridBorders[],
            const IndexType nPoints,
            const int stencilNodes[],
            const ValueType stencilVal[],
            const int stencilOffset[] );
    
        static const char* getId()
        {
            return "Stencil.haloCSR";
        }
    };

    template<typename ValueType>
    struct stencilGEMV
    {
        /** function that computes result += alpha * A * x where A is a stencil matrix.
         *
         *  @param[in,out] result contains values for updated grid points
         *  @param[in] alpha is an additional scaling factor
         *  @param[in] x is the vector which one element for each grid point
         *  @param[in] nDims specifies the dimension of the grid
         *  @param[in] gridSizes contains the dimensions of the grid 
         *  @param[in] width contains maximal left/rights distance for each dimension given by stencil
         *  @param[in] gridDistances contains the distance between two neighbored points for each dim
         *  @param[in] nPoints number of stencil points
         *  @param[in] stencilNodes contins nDims * nPoints direction values for the stencil points
         *  @param[in] stencilVal contains the scale value for each stencil point
         *  @param[in] stencilOffset array with offset for each stencil point in linearized grid
         *
         *  result and x have one value for each grid point corresponding the linearized grid with nDims dimensions.
         *  width might be recomputed by stencilNodes but are passed to avoid a recomputation.
         *  gridDistances might also be recomputed by gridSizes but is used here to be independent from 
         *  column or row-major linearization of grids.
         */
        typedef void ( *FuncType )(
            ValueType result[],
            const ValueType alpha,
            const ValueType x[],
            const IndexType nDims,
            const IndexType gridSizes[],
            const IndexType width[],
            const IndexType gridDistances[],
            const common::Grid::BorderType gridBorders[],
            const IndexType nPoints,
            const int stencilNodes[],
            const ValueType stencilVal[],
            const int stencilOffset[] );

        static const char* getId()
        {
            return "Stencil.GEMV";
        }
    };
};

}

}
