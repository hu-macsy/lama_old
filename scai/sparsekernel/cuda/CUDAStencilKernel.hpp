/**
 * @file CUDAStencilKernel.hpp
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
 * @brief Definition of CUDA routines that implement the stencil kernels.
 * @author Pascal Maczey, Thomas Brandes
 * @date 04.05.2017
 */

#pragma once

#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{
namespace sparsekernel
{

/** Class that contains the CUDA implementations of the stencil kernels. */

class CUDAStencilKernel
{
public:

    /** CUDA implementation for StencilKernelTrait::stencilGEMV */

    template<typename ValueType>
    static void stencilGEMV(
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

private:

    /** Help routine for stencilGEMV on one-dimensional grid */

    template<typename ValueType>
    static void stencilGEMV1(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV2(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV3(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

    template<typename ValueType>
    static void stencilGEMV4(
        ValueType result[],
        const ValueType alpha,
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType nPoints,
        const ValueType stencilVal[],
        const int stencilOffset[] );

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

    CUDAStencilKernel();

    /** Destructor for unregistration. */

    ~CUDAStencilKernel();

    /** Guard for registration during static initialization. */

    static CUDAStencilKernel guard;

    /** Logger for this class. */

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} // namespace sparsekernel

} // namespace scai
