/**
 * @file CUDABLAS1.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief CUDABLAS1.hpp
 * @author lschubert
 * @date 05.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_CUDA_BLAS1_HPP_
#define LAMA_CUDA_BLAS1_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <logging/logging.hpp>

#include <cublas.h>
#include <cuda_runtime_api.h>

namespace lama
{

class SyncToken;    // forward declaration

class LAMA_DLL_IMPORTEXPORT CUDABLAS1
{

// Higher level CUDA BLAS routines of CUDA are allowed to use the private routines
// but routines must not be used directly in LAMA without using the interface

friend class CUDALAPACK;

public:

    /** Routine that sets functions pointers belonging to BLAS1 in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    /**
     * This function is the CUDA implementation of lama::BLASInterface::scal
     */
    template<typename T>
    static void scal( const IndexType n, const T alpha, T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::nrm2
     */
    template<typename T>
    static T nrm2( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::asum
     */
    template<typename T>
    static T asum( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::iamax
     */
    template<typename T>
    static IndexType iamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::swap
     */
    template<typename T>
    static void swap( const IndexType n, T* y, const IndexType incY, T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::copy
     */
    template<typename T>
    static void copy(
        const IndexType n,
        const T* x,
        const IndexType incX,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::axpy
     */
    template<typename T>
    static void axpy(
        const IndexType n,
        const T alpha,
        const T* x,
        const IndexType incX,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::dot
     */
    template<typename T>
    static T dot(
        const IndexType n,
        const T* x,
        const IndexType incX,
        const T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLASInterface::sum
     */
    template<typename T>
    static void sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    template<typename T>
    static void sum_launcher( const int n, T alpha, const T* x, T beta, const T* y, T* z, cudaStream_t stream );

    static bool initialized;   //!< static initialization used for registration

    static bool registerInterface();  //!< registration
};

} /* namespace lama */

#endif // LAMA_CUDA_BLAS1_HPP_
