/**
 * @file CUDABLAS2.hpp
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
 * @brief CUDABLAS2.hpp
 * @author lschubert
 * @date 05.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_CUDABLAS2_HPP_
#define LAMA_CUDABLAS2_HPP_

// for dll_import
#include <scai/common/config.hpp>

// others
#include <scai/lama/LAMATypes.hpp>
#include <scai/tasking/SyncToken.hpp>

// CBLAS_ORDER, CBLAS_TRANSPOSE, ...
#include <scai/lama/cblas.hpp>

// logging
#include <scai/logging.hpp>

#include <cublas_v2.h>
#include <cuda_runtime_api.h>

namespace lama
{

/** Static class that provides CUDA implementaions for the BLAS2 routines of the BLAS interface.
 *
 *  The BLAS2 routines are all private and can only be accessed via registration at an interface.
 *
 */

class COMMON_DLL_IMPORTEXPORT CUDABLAS2
{
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
     * This function is the CUDA implementation of lama::BLAS2Interface::gemv
     */
    template<typename ValueType>
    static void gemv(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE trans,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::symv
     */
    template<typename ValueType>
    static void symv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::trmv
     */
    template<typename ValueType>
    static void trmv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const ValueType* A,
        const IndexType lda,
        ValueType* x,
        const IndexType incX,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::trsv
     */
    template<typename ValueType>
    static void trsv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const ValueType* A,
        const IndexType lda,
        ValueType* x,
        const IndexType incX,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::gbmv
     */
    template<typename ValueType>
    static void gbmv(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE trans,
        const IndexType m,
        const IndexType n,
        const IndexType kl,
        const IndexType ku,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::sbmv
     */
    template<typename ValueType>
    static void sbmv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const IndexType k,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::tbmv
     */
    template<typename ValueType>
    static void tbmv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const IndexType k,
        const ValueType* A,
        const IndexType lda,
        ValueType* x,
        const IndexType incX,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::tbsv
     */
    template<typename ValueType>
    static void tbsv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const IndexType k,
        const ValueType* A,
        const IndexType lda,
        ValueType* x,
        const IndexType incX,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::ger
     */
    template<typename ValueType>
    static void ger(
        const CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        const ValueType* y,
        const IndexType incY,
        ValueType* A,
        const IndexType lda,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::syr
     */
    template<typename ValueType>
    static void syr(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        ValueType* A,
        const IndexType lda,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::syr2
     */
    template<typename ValueType>
    static void syr2(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        const ValueType* y,
        const IndexType incY,
        ValueType* A,
        const IndexType lda,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::spr
     */
    template<typename ValueType>
    static void spr(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        ValueType* AP,
        tasking::SyncToken* syncToken );

    /**
     * This function is the CUDA implementation of lama::BLAS2Interface::spr2
     */
    template<typename ValueType>
    static void spr2(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        const ValueType* y,
        const IndexType incY,
        ValueType* AP,
        tasking::SyncToken* syncToken );

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static    bool initialized; //!< static initialization used for registration

    static bool registerInterface();//!< registration
};

}
/* namespace lama */

#endif // LAMA_CUDABLAS2_HPP_
