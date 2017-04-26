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

class OpenMPStencilKernel
{
public:

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
        const int stencilGridPos[], 
        const ValueType stencilVal[],
        const int stencilLinPos[] );

private:

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
        const int stencilGridPos[], 
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
