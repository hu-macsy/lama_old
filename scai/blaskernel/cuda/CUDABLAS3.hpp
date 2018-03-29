/**
 * @file CUDABLAS3.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of static class with CUDA implementations of BLAS3 routines.
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// local library
#include <scai/blaskernel/cblas.hpp>

// internal scai library
#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

#include <scai/logging.hpp>

#include <scai/kregistry/mepr/Registrator.hpp>

namespace scai
{

namespace blaskernel
{

/** Static class that provides CUDA implementaions for the BLAS3 routines as specified in BLASKernelTrait.
 *
 *  The BLAS3 routines are all private and can only be accessed via kernel registry.
 */

class COMMON_DLL_IMPORTEXPORT CUDABLAS3
{
private:

    /**
     * @brief CUDA implementation of BLASKernelTrait::gemm
     * 
     * This implementation will call a corresponding cuBLAS function.
     */
    template<typename ValueType>
    static void gemm(
        const common::MatrixOp opA,
        const common::MatrixOp opB,
        const IndexType m,
        const IndexType n,
        const IndexType k,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        const ValueType* B,
        const IndexType ldb,
        const ValueType beta,
        ValueType* C,
        const IndexType ldc );

    /**
     * @brief trsm solves one of the matrix equations
     *
     * op(A) * X = alpha * B or
     * X * op(A) = alpha * B
     * where op(A) = A or op(A) = AT
     *
     * alpha is a scalar,
     * X and B are m×n matrices that consist of elements.
     * A is a unit or non‐unit, upper or lower, triangular matrix.
     * The result matrix X overwrites input matrix B;
     * that is, on exit the result is stored in B.
     * matrix A and B are stored in column‐major format,
     * lda and ldb are the leading dimensions of the two‐dimensional arrays
     * that contain A and B, respectively.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] side    specifies whether op(A) appears on the left or right of X:
     *                    If side == 'L' or 'l', op(A) * X = alpha * B
     *                    If side == 'R' or 'R', X * op(A) = alpha * B
     * @param[in] uplo    specifies whether the  matrix A is an upper or lower triangular matrix.
     *                    If uplo == 'U' or 'u',
     *                    A is an upper triangular matrix.
     *                    If uplo == 'L' or 'l',
     *                    A is a lower triangular matrix.
     * @param[in] transa  specifies the form of op(A) to be used in the matrix multiplication.
     *                    If transa == 'N' or 'n', op(A) = A.
     *                    If transa == 'T','t','C','c', op(A) = AT.
     * @param[in] diag    specifies whether or not A is a unit triangular matrix.
     *                    If diag == 'U' or 'u',
     *                    A is assumed to be unit triangular.
     *                    If diag == 'N' or 'n',
     *                    A is not assumed to be unit triangular.
     * @param[in] m       the number of rows of matrix B;
     *                    m must be at least zero.
     * @param[in] n       the number of columns of matrix B;
     *                    n must be at least zero.
     * @param[in] alpha   scalar multiplier applied to B.
     *                    When alpha is zero,
     *                    A is not referenced and B does not have to be a valid input.
     * @param[in] A       array of dimensions (lda, k),
     *                    where k is m when side == 'L' or 'l'
     *                    and is n when side == 'R' or 'r'.
     *                    If uplo == 'U' or 'u',
     *                    the leading k×k upper triangular part of the array A
     *                    must contain the upper triangular matrix,
     *                    and the strictly lower triangular matrix of A is not referenced.
     *                    When uplo == 'L' or 'l',
     *                    the leading k×k lower triangular part of the array A
     *                    must contain the lower triangular matrix,
     *                    and the strictly upper triangular part of A is not referenced.
     *                    If diag == 'U' or 'u',
     *                    the diagonal elements of A are not referenced and are assumed to be unity.
     * @param[in] lda     leading dimension of the two-dimensional array containing A.
     *                    When side == 'L' or 'l', lda must be at least max(1, m).
     *                    When side == 'R' or 'r', lda must be at least max(1, n).
     * @param[in] B       array of dimensions (ldb, n);
     *                    ldb must be at least max(1, m). ValueType
     *                    he leading m×n part of the array B must contain
     *                    the righthand side matrix B.
     *                    On exit B is overwritten by the solution matrix X.
     * @param[in] ldb     leading dimension of the two-dimensional array containing B;
     *                    ldb must be at least max(1, m).
     */
    template<typename ValueType>
    static void trsm(
        const CBLAS_ORDER order,
        const CBLAS_SIDE side,
        const CBLAS_UPLO uplo,
        const CBLAS_TRANSPOSE transa,
        const CBLAS_DIAG diag,
        const IndexType m,
        const IndexType n,
        const ValueType alpha,
        const ValueType* A,
        const IndexType lda,
        ValueType* B,
        const IndexType ldb );

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    /** Registration of methods  at kernel registry. */

    template<typename ValueType>
    struct RegistratorV
    {
        static void registerKernels( const kregistry::KernelRegistry::KernelRegistryFlag flag );
    };

    /** Constructor for registration. */

    CUDABLAS3();

    /** Destructor for unregistration. */

    ~CUDABLAS3();

    /** Static variable for registration at static initialization. */

    static CUDABLAS3 guard;

};

} /* end namespace blaskernel */

} /* end namespace scai */
