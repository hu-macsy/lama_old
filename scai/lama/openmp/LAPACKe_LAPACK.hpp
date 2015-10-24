/**
 * @file LAPACKe_LAPACK.hpp
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
 * @brief Implemenation of LAMA LAPACK interface via LAPACKe
 * @author Thomas Brandes
 * @date 10.04.2013
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/openmp/OpenMP.hpp>
#include <scai/lama/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

namespace scai
{

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

    /** Implementation of BLASKernelTrait::LAPACK::getrf by LAPACKe. */

    template<typename ValueType>
    static IndexType getrf(
        const CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        ValueType* const a,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getri by LAPACKe. */

    template<typename ValueType>
    static IndexType getri(
        const CBLAS_ORDER order,
        const IndexType n,
        ValueType* const A,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getinv by LAPACKe. */

    template<typename ValueType>
    static void getinv( const IndexType n, ValueType* a, const IndexType lda );

    /** Implementation of BLASKernelTrait::LAPACK::tptrs by LAPACKe. */

    template<typename ValueType>
    static IndexType tptrs(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const IndexType nrhs,
        const ValueType* AP,
        ValueType* B,
        const IndexType ldb );

    /** Routine that sets functions pointers belonging to LAPACK in a BLASKernelTrait.
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void registerKernels();

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    static    bool initialized;

    static bool registerInterface();

}; /* LAPACKe_LAPACK */

} /* end namespace lama */

} /* end namespace scai */
