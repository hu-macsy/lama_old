/**
 * @file CUSolverCSRUtils.hpp
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
 * @brief Provide CSR routines by using CUSolver library
 * @author Lauretta Schubert
 * @date 19.07.2016
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/kregistry/mepr/Registrator.hpp>

#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/assert.hpp>

#include <cuda_runtime_api.h>

#ifndef CUDART_VERSION
#error CUDART_VERSION Undefined!
#elif ( CUDART_VERSION >= 7050 )
#include <cusolverSp.h>

namespace scai
{

namespace sparsekernel
{

class COMMON_DLL_IMPORTEXPORT CUSolverCSRUtils
{
public:

    /** Implementation for CSRKernelTrait::decomposition */

    template<typename ValueType>
    static void decomposition(
        ValueType* const solution,
        const IndexType csrIA[],
        const IndexType csrJA[],
        const ValueType csrValues[],
        const ValueType rhs[],
        const IndexType numRows,
        const IndexType nnz,
        const bool isSymmetic );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

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

    CUSolverCSRUtils();

    /** Destructor for unregistration. */

    ~CUSolverCSRUtils();

    /** Static variable for registration at static initialization. */

    static CUSolverCSRUtils guard;
};

/* --------------------------------------------------------------------------- */

} /* end namespace sparsekernel */

} /* end namespace scai */

#endif
