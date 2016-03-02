/**
 * @file LAPACK_LAPACK.hpp
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
 * @brief LAPACK_LAPACK.hpp
 * @author Lauretta Schubert
 * @date 02.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/tasking/SyncToken.hpp>

#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/kregistry/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

class LAPACK_LAPACK
{
public:

    /** Implementation of BLASKernelTrait::LAPACK::getrf by LAPACK. */

    template<typename ValueType>
    static IndexType getrf(
        const CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        ValueType* const a,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getri by LAPACK. */

    template<typename ValueType>
    static IndexType getri(
        const CBLAS_ORDER order,
        const IndexType n,
        ValueType* const A,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of BLASKernelTrait::LAPACK::getinv by LAPACK. */

    template<typename ValueType>
    static void getinv( const IndexType n, ValueType* a, const IndexType lda );

    /** Implementation of BLASKernelTrait::LAPACK::tptrs vi LAPACK. */

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

    /** OpenMP implementation for BLASKernelTrait::LAPACK::laswp */

    template<typename ValueType>
    static void laswp(
        const CBLAS_ORDER order,
        const IndexType n,
        ValueType* A,
        const IndexType lda,
        const IndexType k1,
        const IndexType k2,
        const IndexType* ipiv,
        const IndexType incx );

private:


    /** Routine that registers all methods at the kernel registry. */

    SCAI_DECLARE_REGISTRATOR( Registrator, typename ValueType )

    /** Constructor for registration. */

    LAPACK_LAPACK();

    /** Destructor for unregistration. */

    ~LAPACK_LAPACK();

    /** Static variable for registration at static initialization. */

    static LAPACK_LAPACK guard;

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

}; /* LAPACK_LAPACK */

} /* end namespace blaskernel */

} /* end namespace scai */
