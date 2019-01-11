/**
 * @file CUDABLAS1.hpp
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
 * @brief Definition of static class with CUDA implementations of BLAS1 routines.
 * @author Eric Schricker
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

// CUDA
#include <cublas_v2.h>
#include <cuda_runtime_api.h>

namespace scai
{

namespace blaskernel
{

/** Static class that provides CUDA implementaions for the BLAS1 routines as specified in BLASKernelTrait.
 *
 *  The BLAS1 routines are all private and can only be accessed via kernel registry.
 */

class COMMON_DLL_IMPORTEXPORT CUDABLAS1
{

// Higher level CUDA BLAS routines of CUDA are allowed to use the private routines
// but routines must not be used directly in LAMA without using the interface

    friend class CUDALAPACK;

private:

    /**
     * This function is the CUDA implementation of BLASKernelTrait::scal
     */
    template<typename ValueType>
    static void scal(
        const IndexType n,
        const ValueType alpha,
        ValueType* x,
        const IndexType incX );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::nrm2
     */
    template<typename ValueType>
    static ValueType nrm2( const IndexType n, const ValueType* x, const IndexType incX );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::asum
     */
    template<typename ValueType>
    static ValueType asum( const IndexType n, const ValueType* x, const IndexType incX );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::iamax
     */
    template<typename ValueType>
    static IndexType iamax( const IndexType n, const ValueType* x, const IndexType incX );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::swap
     */
    template<typename ValueType>
    static void swap(
        const IndexType n,
        ValueType* y,
        const IndexType incY,
        ValueType* x,
        const IndexType incX );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::copy
     */
    template<typename ValueType>
    static void copy(
        const IndexType n,
        const ValueType* x,
        const IndexType incX,
        ValueType* y,
        const IndexType incY );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::axpy
     */
    template<typename ValueType>
    static void axpy(
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        ValueType* y,
        const IndexType incY );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::dot
     */
    template<typename ValueType>
    static ValueType dot(
        const IndexType n,
        const ValueType* x,
        const IndexType incX,
        const ValueType* y,
        const IndexType incY );

    /**
     * This function is the CUDA implementation of BLASKernelTrait::sum
     */
    template<typename ValueType>
    static void sum(
        const IndexType n,
        ValueType alpha,
        const ValueType* x,
        ValueType beta,
        const ValueType* y,
        ValueType* z );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    template    <typename ValueType>
    static void sum_launcher( const int n, ValueType alpha, const ValueType* x, ValueType beta, const ValueType* y, ValueType* z, cudaStream_t stream );

    /** Register or unregister all kernel implementations of this class. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDABLAS1();

    /** Destructor for unregistration. */

    ~CUDABLAS1();

    /** Static variable for registration at static initialization. */

    static CUDABLAS1 guard;
};

} /* end namespace blaskernel */

} /* end namespace scai */
