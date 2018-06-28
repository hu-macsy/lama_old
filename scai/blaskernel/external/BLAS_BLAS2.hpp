/**
 * @file BLAS_BLAS2.hpp
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
 * @brief BLAS2 utilities on Host Context by wrapping to BLAS library
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** Implementations of blas2 methods of BLASKernelTrait with OpenMP.  */

class COMMON_DLL_IMPORTEXPORT BLAS_BLAS2
{
public:

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::gemv
     */
    template<typename ValueType>
    static void gemv(
        const common::MatrixOp opA,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY );

private:

    /** Routine that registers all methods at the kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    BLAS_BLAS2();

    /** Destructor for unregistration. */

    ~BLAS_BLAS2();

    /** Static variable for registration at static initialization. */

    static BLAS_BLAS2 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */
