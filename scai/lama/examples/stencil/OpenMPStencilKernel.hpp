/**
 * @file OpenMPStencilKernel.hpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Definition of OpenMP routines that implement the stencil kernels.
 * @author Pascal Maczey, Thomas Brandes
 * @date Apr 26, 2017
 */

#pragma once

#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{
namespace stencilkernel
{

/** Class that contains the OpenMP implementations of the stencil kernels. */

class OpenMPStencilKernel
{
public:

    /** OpenMP implementation for StencilKernelTrait::stencilSizes */

    static void stencilSizes(
        IndexType sizes[],
        const IndexType nDims, 
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[] );

    /** OpenMP implementation for StencilKernelTrait::stencil2CSR */

    template<typename ValueType>
    static void stencil2CSR(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType nDims,
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[],
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** OpenMP implementation for StencilKernelTrait::stencilGEMV */

    template<typename ValueType>
    static void stencilGEMV(
        ValueType result[], 
        const ValueType alpha,  
        const ValueType x[],
        const IndexType nDims, 
        const IndexType gridSizes[],
        const IndexType lb[],
        const IndexType ub[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[], 
        const ValueType stencilVal[],
        const int stencilLinPos[] );

private:

    /** Implementation of stencilSizes for nDims == 1 */

    static void stencilSizes1(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[] );

    /** Implementation of stencilSizes for nDims == 2 */

    static void stencilSizes2(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[] );

    /** Implementation of stencilSizes for nDims == 3 */

    static void stencilSizes3(
        IndexType sizes[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[] );

    /** Implementation of stencil2CSR for nDims == 1 */

    template<typename ValueType>
    static void stencil2CSR1(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[],
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** Implementation of stencil2CSR for nDims == 2 */

    template<typename ValueType>
    static void stencil2CSR2(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[],
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** Implementation of stencil2CSR for nDims == 3 */

    template<typename ValueType>
    static void stencil2CSR3(
        IndexType csrJA[],
        ValueType csrValues[],
        const IndexType csrIA[],
        const IndexType gridSizes[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[],
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** Implementation of stencilGEMV for nDims == 1 */

    template<typename ValueType>
    static void stencilGEMV1(
        ValueType result[], 
        const ValueType alpha,  
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType lb[],
        const IndexType ub[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[], 
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** Implementation of stencilGEMV for nDims == 2 */

    template<typename ValueType>
    static void stencilGEMV2(
        ValueType result[], 
        const ValueType alpha,  
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType lb[],
        const IndexType ub[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[], 
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** Implementation of stencilGEMV for nDims == 3 */

    template<typename ValueType>
    static void stencilGEMV3(
        ValueType result[], 
        const ValueType alpha,  
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType lb[],
        const IndexType ub[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[], 
        const ValueType stencilVal[],
        const int stencilLinPos[] );

    /** Implementation of stencilGEMV for nDims == 4 */

    template<typename ValueType>
    static void stencilGEMV4(
        ValueType result[], 
        const ValueType alpha,  
        const ValueType x[],
        const IndexType gridSizes[],
        const IndexType lb[],
        const IndexType ub[],
        const IndexType gridDistances[],
        const IndexType nPoints,
        const int stencilNodes[], 
        const ValueType stencilVal[],
        const int stencilLinPos[] );

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

} // namespace stencilkernel

} // namespace scai
