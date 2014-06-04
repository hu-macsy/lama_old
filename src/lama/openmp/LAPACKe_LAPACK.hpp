/**
 * @file LAPACKe_LAPACK.hpp
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
 * @brief Implemenation of LAMA LAPACK interface via LAPACKe
 * @author Thomas Brandes
 * @date 10.04.2013
 * @since 1.0.0
 */
#ifndef LAMA_OPENMP_LAPACKE_HPP_
#define LAMA_OPENMP_LAPACKE_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>
#include <lama/openmp/OpenMP.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** @brief Implementation of LAPACK functionality at the LAMAInterface by
 *         using LAPACKe (C/C++ extensions for LAPACK).
 *
 *  The use of these routines is more convenient if LAPACKe is available
 *  but has the same functionality as OpenMPLAPACK using Fortran calling conventions.
 */

class LAPACKe_LAPACK
{
public:

    /** Implementation of BLASInterface::LAPACK::getrf by LAPACKe. */

    template<typename T>
    static IndexType getrf(
        const enum CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        T* const a,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASInterface::LAPACK::getri by LAPACKe. */

    template<typename T>
    static IndexType getri(
        const enum CBLAS_ORDER order,
        const IndexType n,
        T* const A,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASInterface::LAPACK::getinv by LAPACKe. */

    template<typename T>
    static void getinv( const IndexType n, T* a, const IndexType lda );

    /** Implementation of BLASInterface::LAPACK::tptrs by LAPACKe. */

    template<typename T>
    static IndexType tptrs(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const IndexType nrhs,
        const T* AP,
        T* B,
        const IndexType ldb );

    /** Routine that sets functions pointers belonging to LAPACK in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in OpenMP
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static bool initialized;

    static bool registerInterface();

}; /* LAPACKe_LAPACK */

} /* namespace lama */

#endif // LAMA_OPENMP_LAPACKE_HPP_
