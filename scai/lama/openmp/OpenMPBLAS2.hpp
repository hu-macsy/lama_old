/**
 * @file OpenMPBLAS2.hpp
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
 * @brief OpenMPBLAS2.hpp
 * @author Eric Schricker
 * @date 05.07.2012
 * @since 1.0.0
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/lama/openmp/BLASHelper.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <scai/common/SCAITypes.hpp>

namespace scai
{

namespace lama
{

/** Implementations of methods for scai::lama::BLAS2Interface with OpenMP.
 *
 *  @todo Move all method documentations to LAMAInterface and make references here
 *  @todo Add information here about use of native BLAS2 libraries
 */

class COMMON_DLL_IMPORTEXPORT OpenMPBLAS2
{
public:

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::gemv
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
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::trmv
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
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::trsv
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
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::gbmv
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
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::sbmv
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
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::tbmv
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
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::tbsv
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
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::ger
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
        const IndexType lda );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::syr
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
        const IndexType lda );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::syr2
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
        const IndexType lda );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::spmv
     */
    template<typename ValueType>
    static void spmv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* AP,
        const ValueType* x,
        const IndexType incX,
        const ValueType beta,
        ValueType* y,
        const IndexType incY );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::spr
     */
    template<typename ValueType>
    static void spr(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const IndexType n,
        const ValueType alpha,
        const ValueType* x,
        const IndexType incX,
        ValueType* AP );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::spr2
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
        ValueType* AP );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::tpmv
     */
    template<typename ValueType>
    static void tpmv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const ValueType* AP,
        ValueType* x,
        const IndexType incX );

    /**
     * This function is the OpenMP implementation of scai::lama::BLAS2Interface::tpsv
     */
    template<typename ValueType>
    static void tpsv(
        const CBLAS_ORDER order,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE trans,
        const CBLAS_DIAG diag,
        const IndexType n,
        const ValueType* Ap,
        ValueType* x,
        const IndexType incX );

    /**
     * @todo add doxygen comment
     * @todo clarify BLAS inteface
     */
    template<typename ValueType>
    static void agemvpbv(
        int n,
        const double alpha,
        const double* const a,
        int m,
        const double* const x,
        const double beta,
        const double* const z,
        double* y );

    /** Routine that sets functions pointers belonging to BLAS2 in a BLASKernelTrait.
     *
     *  param[inout] BLASKernelTrait struct to register all routines implemented in CUDA
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void registerKernels();

private:

    static bool initialized;

    static bool registerInterface();

    SCAI_LOG_DECL_STATIC_LOGGER( logger )
};

} /* end namespace lama */

} /* end namespace scai */
