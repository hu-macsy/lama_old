/**
 * @file OpenMPBLAS1.hpp
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
 * @brief OpenMPBLAS1.hpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_OPENMPBLAS1_HPP_
#define LAMA_OPENMPBLAS1_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** Implementations of methods for lama::BLAS1Interface with OpenMP.
 *
 *  @todo Add information here about use of native BLAS1 libraries
 */

class LAMA_DLL_IMPORTEXPORT OpenMPBLAS1
{
public:

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::BLAS1::scal
     */
    template<typename T>
    static void scal( const IndexType n, const T alpha, T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::nrm2
     */
    template<typename T>
    static T nrm2( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::asum
     */
    template<typename T>
    static T asum( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::iamax
     */
    template<typename T>
    static IndexType iamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::viamax
     */
    template<typename T>
    static T viamax( const IndexType n, const T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::swap
     */
    template<typename T>
    static void swap( const IndexType n, T* y, const IndexType incY, T* x, const IndexType incX, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::copy
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
     * This function is the OpenMP implementation of lama::BLAS1Interface::axpy
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
     * This function is the OpenMP implementation of lama::BLAS1Interface::dot
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
     * This function is the OpenMP implementation of lama::BLAS1Interface::sum
     */
    template<typename T>
    static void sum( const IndexType n, T alpha, const T* x, T beta, const T* y, T* z, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::rot
     */
    template<typename T>
    static void rot(
        const IndexType n,
        T* x,
        const IndexType incX,
        T* y,
        const IndexType incY,
        const T c,
        const T s,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::rotm
     */
    template<typename T>
    static void rotm(
        const IndexType n,
        T* x,
        const IndexType incX,
        T* y,
        const IndexType incY,
        const T* P,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::rotg
     */
    template<typename T>
    static void rotg( T* a, T* b, T* c, T* s, SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS1Interface::rotmg
     */
    template<typename T>
    static void rotmg( T* d1, T* d2, T* b1, const T b2, T* P, SyncToken* syncToken );

    /** Routine that sets functions pointers belonging to BLAS1 in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    static bool initialized;

    static bool registerInterface();

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

} /* namespace lama */

#endif // LAMA_OPENMPBLAS1_HPP_
