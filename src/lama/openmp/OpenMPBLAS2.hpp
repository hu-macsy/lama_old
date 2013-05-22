/**
 * @file OpenMPBLAS2.hpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief OpenMPBLAS2.hpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * $Id$
 */
#ifndef LAMA_OPENMPBLAS2_HPP_
#define LAMA_OPENMPBLAS2_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>

#include <lama/openmp/BLASHelper.hpp>

// logging
#include <logging/logging.hpp>

namespace lama
{

/** Implementations of methods for lama::BLAS2Interface with OpenMP.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS1 libraries
 */

class LAMA_DLL_IMPORTEXPORT OpenMPBLAS2
{
public:

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::gemv
     */
    template<typename T>
    static void gemv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_TRANSPOSE trans,
        const IndexType m,
        const IndexType n,
        const T alpha,
        const T* A,
        const IndexType lda,
        const T* x,
        const IndexType incX,
        const T beta,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::symv
     */
    template<typename T>
    static void symv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const T alpha,
        const T* A,
        const IndexType lda,
        const T* x,
        const IndexType incX,
        const T beta,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::trmv
     */
    template<typename T>
    static void trmv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const T* A,
        const IndexType lda,
        T* x,
        const IndexType incX,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::trsv
     */
    template<typename T>
    static void trsv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const T* A,
        const IndexType lda,
        T* x,
        const IndexType incX,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::gbmv
     */
    template<typename T>
    static void gbmv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_TRANSPOSE trans,
        const IndexType m,
        const IndexType n,
        const IndexType kl,
        const IndexType ku,
        const T alpha,
        const T* A,
        const IndexType lda,
        const T* x,
        const IndexType incX,
        const T beta,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::sbmv
     */
    template<typename T>
    static void sbmv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const IndexType k,
        const T alpha,
        const T* A,
        const IndexType lda,
        const T* x,
        const IndexType incX,
        const T beta,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::tbmv
     */
    template<typename T>
    static void tbmv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const IndexType k,
        const T* A,
        const IndexType lda,
        T* x,
        const IndexType incX,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::tbsv
     */
    template<typename T>
    static void tbsv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const IndexType k,
        const T* A,
        const IndexType lda,
        T* x,
        const IndexType incX,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::ger
     */
    template<typename T>
    static void ger(
        const enum CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        const T alpha,
        const T* x,
        const IndexType incX,
        const T* y,
        const IndexType incY,
        T* A,
        const IndexType lda,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::syr
     */
    template<typename T>
    static void syr(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const T alpha,
        const T* x,
        const IndexType incX,
        T* A,
        const IndexType lda,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::syr2
     */
    template<typename T>
    static void syr2(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const T alpha,
        const T* x,
        const IndexType incX,
        const T* y,
        const IndexType incY,
        T* A,
        const IndexType lda,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::spmv
     */
    template<typename T>
    static void spmv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const T alpha,
        const T* AP,
        const T* x,
        const IndexType incX,
        const T beta,
        T* y,
        const IndexType incY,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::spr
     */
    template<typename T>
    static void spr(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const T alpha,
        const T* x,
        const IndexType incX,
        T* AP,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::spr2
     */
    template<typename T>
    static void spr2(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const IndexType n,
        const T alpha,
        const T* x,
        const IndexType incX,
        const T* y,
        const IndexType incY,
        T* AP,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::tpmv
     */
    template<typename T>
    static void tpmv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const T* AP,
        T* x,
        const IndexType incX,
        SyncToken* syncToken );

    /**
     * This function is the OpenMP implementation of lama::BLAS2Interface::tpsv
     */
    template<typename T>
    static void tpsv(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const T* Ap,
        T* x,
        const IndexType incX,
        SyncToken* syncToken );

    /**
     * @todo add doxygen comment
     * @todo clarify BLAS inteface
     */
    template<typename T>
    static void agemvpbv(
        int n,
        const double alpha,
        const double* const a,
        int m,
        const double* const x,
        const double beta,
        const double* const z,
        double* y,
        SyncToken* syncToken );

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

#endif // LAMA_OPENMPBLAS2_HPP_
