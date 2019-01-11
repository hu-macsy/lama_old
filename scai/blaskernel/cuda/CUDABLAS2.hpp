/**
 * @file CUDABLAS2.hpp
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
 * @brief Definition of static class with CUDA implementations of BLAS2 routines.
 * @author lschubert
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp> // CBLAS_ORDER, CBLAS_TRANSPOSE, ...

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

// CUDA
#include <cublas_v2.h>
#include <cuda_runtime_api.h>

namespace scai
{

namespace blaskernel
{

/** Static class that provides CUDA implementaions for the BLAS2 routines as specified in BLASKernelTrait.
 *
 */

class COMMON_DLL_IMPORTEXPORT CUDABLAS2
{
public:

    /**
     * This function is the CUDA implementation of BLASKernelTrait::gemv
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

    /**
     * This function is the CUDA implementation of BLASKernelTrait::geam
     */
    template<typename ValueType>
    static void geam(
        ValueType* C,
        const IndexType ldc,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const common::MatrixOp opA,
        const ValueType beta,
        const ValueType* B,
        const IndexType ldb,
        const common::MatrixOp opB );

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** Registration of methods at kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDABLAS2();

    /** Destructor for unregistration. */

    ~CUDABLAS2();

    /** Static variable for registration at static initialization. */

    static CUDABLAS2 guard;

};

} /* end namespace blaskernel */

} /* end namespace scai */
