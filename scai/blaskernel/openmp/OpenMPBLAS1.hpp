/**
 * @file OpenMPBLAS1.hpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Class with implementation of BLAS1 routines in C++ with OpenMP parallelization.
 * @author Eric Schricker
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** Implementations of methods for BLASKernelTrait with OpenMP.
 *
 *  Instead of using native BLAS1 libraries this class has own C++
 *  implementation that are portable across all platforms.
 *  Furthermore, due to OpenMP parallelization, these routines are
 *  rather fast.
 */

class COMMON_DLL_IMPORTEXPORT OpenMPBLAS1
{
public:

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::scal
     */
    template<typename ValueType>
    static void scal(
        const IndexType n,
        const ValueType alpha,
        ValueType* x,
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::nrm2
     */
    template<typename ValueType>
    static ValueType nrm2( const IndexType n, const ValueType* x, const IndexType incX );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::asum
     */
    template<typename ValueType>
    static ValueType asum( const IndexType n, const ValueType* x, const IndexType incX );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::iamax
     */
    template<typename ValueType>
    static IndexType iamax( const IndexType n, const ValueType* x, const IndexType incX );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::swap
     */
    template<typename ValueType>
    static void swap(
        const IndexType n,
        ValueType* y,
        const IndexType incY,
        ValueType* x,
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::copy
     */
    template<typename ValueType>
    static void copy(
        const IndexType n,
        const ValueType* x,
        const IndexType incX,
        ValueType* y,
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::axpy
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
     * This function is the OpenMP implementation of BLASKernelTrait::dot
     */
    template<typename ValueType>
    static ValueType dot(
        const IndexType n,
        const ValueType* x,
        const IndexType incX,
        const ValueType* y,
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of BLASKernelTrait::sum
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

    /** structure that registers all methods at the kernel registry. */

    SCAI_DECLARE_REGISTRATOR( RegistratorV, template<typename ValueType> )

    /** Constructor for registration. */

    OpenMPBLAS1();

    /** Destructor for unregistration. */

    ~OpenMPBLAS1();

    /** Static variable for registration at static initialization. */

    static OpenMPBLAS1 guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace blaskernel */

} /* end namespace scai */
