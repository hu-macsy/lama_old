/**
 * @file CUSolverCSRUtils.hpp
 *
 * @license
 * Copyright (c) 2009-2016
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
        const IndexType nnz );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Routine that registers all methods at the kernel registry. */

    SCAI_KREGISTRY_DECL_REGISTRATOR( Registrator )
    SCAI_KREGISTRY_DECL_REGISTRATOR( RegistratorV, template<typename ValueType> )

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
