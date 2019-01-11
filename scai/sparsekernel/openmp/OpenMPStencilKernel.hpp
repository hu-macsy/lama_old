/**
 * @file OpenMPStencilKernel.hpp
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
 * @brief Definition of OpenMP routines that implement the stencil kernels.
 * @author Pascal Maczey, Thomas Brandes
 * @date 26.04.2017
 */

#pragma once

#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>
#include <scai/common/Grid.hpp>

namespace scai
{
namespace sparsekernel
{

/** Class that contains the OpenMP implementations of the stencil kernels. */

class COMMON_DLL_IMPORTEXPORT OpenMPStencilKernel
{
public:

    /** OpenMP implementation for StencilKernelTrait::stencilLocalSizes */

    static void stencilLocalSizes(
        IndexType sizes[],
        const IndexType nDims, 
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** OpenMP implementation for StencilKernelTrait::stencilLocalCSR */

    template<typename ValueType>
    static void stencilLocalCSR(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType nDims,
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** OpenMP implementation for StencilKernelTrait::stencilHaloSizes */

    static void stencilHaloSizes(
        IndexType sizes[],
        const IndexType nDims,
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** OpenMP implementation for StencilKernelTrait::stencilHaloCSR */

    template<typename ValueType>
    static void stencilHaloCSR(
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
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** OpenMP implementation for StencilKernelTrait::normalGEMV */

    template<typename ValueType>
    static void normalGEMV(
        ValueType result[], 
        const ValueType alpha,  
        const ValueType x[],
        const ValueType beta,  
        const ValueType y[],
        const IndexType nDims, 
        const IndexType hostGridSizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType gridStencilWidth[],
        const IndexType nPoints,
        const int stencilPositions[], 
        const ValueType stencilVal[],
        const int stencilOffset[] );

private:

    /** Implementation of stencilLocalSizes for nDims == 1 */

    static void stencilLocalSizes1(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilLocalSizes for nDims == 2 */

    static void stencilLocalSizes2(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilLocalSizes for nDims == 3 */

    static void stencilLocalSizes3(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilLocalSizes for nDims == 4 */

    static void stencilLocalSizes4(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilLocalCSR for nDims == 1 */

    template<typename ValueType>
    static void stencilLocalCSR1(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilLocalCSR for nDims == 2 */

    template<typename ValueType>
    static void stencilLocalCSR2(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilLocalCSR for nDims == 3 */

    template<typename ValueType>
    static void stencilLocalCSR3(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilLocalCSR for nDims == 4 */

    template<typename ValueType>
    static void stencilLocalCSR4(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const common::Grid::BorderType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilHaloSizes for nDims == 1 */

    static void stencilHaloSizes1(
        IndexType sizes[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilHaloSizes for nDims == 2 */

    static void stencilHaloSizes2(
        IndexType sizes[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilHaloSizes for nDims == 3 */

    static void stencilHaloSizes3(
        IndexType sizes[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilHaloSizes for nDims == 4 */

    static void stencilHaloSizes4(
        IndexType sizes[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[] );

    /** Implementation of stencilHaloCSR for nDims == 1 */

    template<typename ValueType>
    static void stencilHaloCSR1(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const IndexType globalGridDistances[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilHaloCSR for nDims == 2 */

    template<typename ValueType>
    static void stencilHaloCSR2(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const IndexType globalGridDistances[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilHaloCSR for nDims == 3 */

    template<typename ValueType>
    static void stencilHaloCSR3(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const IndexType globalGridDistances[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Implementation of stencilHaloCSR for nDims == 4 */

    template<typename ValueType>
    static void stencilHaloCSR4(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType localGridSizes[],
        const IndexType localGridDistances[],
        const IndexType localLB[],
        const IndexType globalGridSizes[],
        const IndexType globalGridDistances[],
        const common::Grid::BorderType globalGridBorders[],
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV1Inner(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV2Inner(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV3Inner(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV4Inner(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMVInner(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType nDims,
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV1Border(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[], 
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV2Border(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[], 
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV3Border(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[], 
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV4Border(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[], 
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMVBorder(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType nDims,
        const IndexType gridSizes[],
        const IndexType gridBounds[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType nPoints,
        const int stencilPositions[], 
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMVCaller(
        IndexType gridBounds[],
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType nDims,
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType gridBorders[],
        const IndexType gridStencilwidth[],
        const IndexType currentDim,
        const IndexType nPoints,
        const int stencilPositions[],
        const ValueType stencilVal[],
        const int stencilOffset[] );

    /** Struct for registration of methods without template arguments */

    struct Registrator
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Struct for registration of methods with one template argument.
     *
     *  Registration function is wrapped in struct/class that can be used as template
     *  argument for metaprogramming classes to expand for each supported type
     */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    OpenMPStencilKernel();

    /** Destructor for unregistration. */

    ~OpenMPStencilKernel();

    /** Guard for registration during static initialization. */

    static OpenMPStencilKernel guard;

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace sparsekernel

} // namespace scai
