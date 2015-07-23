/**
 * @file CUDABLAS3.hpp
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
 * @brief CUDABLAS3.hpp
 * @author lschubert
 * @date 05.07.2012
 * @since 1.0.0
 */
#ifndef LAMA_CUDABLAS3_HPP_
#define LAMA_CUDABLAS3_HPP_

// for dll_import
#include <common/config.hpp>

// others
#include <lama/LAMATypes.hpp>

// logging
#include <logging/logging.hpp>

#include <lama/cblas.hpp>

namespace tasking
{
    class SyncToken;
}

namespace lama
{

class COMMON_DLL_IMPORTEXPORT CUDABLAS3
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
     * @brief gemm computes the product of matrix A and matrix B,
     * multiplies the result by scalar alpha,
     * and adds the sum to the product of matrix C and scalar beta.
     * It performs one of the matrix‐matrix operations:
     *
     * C = alpha * op(A) * op(B) + beta * C,
     * where op(X) = X or op(X) = XT
     *
     * alpha and beta are scalars.
     * A, B, and C are matrices consisting of elements,
     * with op(A) an m×k matrix,
     * op(B) a k×n matrix,
     * C an m×n matrix.
     * matrix A, B, and C are stored in column‐major format,
     * lda, ldb, and ldc are the leading dimensions of the two‐dimensional arrays
     * containing A, B, and C.
     *
     * @param[in] order   TODO[doxy] Complete Description.
     * @param[in] transa  specifies op(A)
     *                    If transa = 'N', 'n'
     *                    op(A) = A
     *                    If transa == 'T', 't', 'C', 'c'
     *                    op(A) = AT
     * @param[in] transb  specifies op(B)
     *                    If transb = 'N', 'n'
     *                    op(B) = B
     *                    If transb == 'T', 't', 'C', 'c'
     *                    op(B) = BT
     * @param[in] m       number of rows of matrix op(A) and rows of matrix C;
     *                    m must be at least zero.
     * @param[in] n       number of columns of matrix op(B) and number of columns of C;
     *                    n must be at least zero.
     * @param[in] k       number of columns of matrix op(A) and number of rows of op(B);
     *                    k must be at least zero.
     * @param[in] alpha   scalar multiplier applied to op(A) * op(B)
     * @param[in] A       array of dimensions (lda, k)
     *                    if transa == 'N' or 'n', and of dimensions (lda, m) otherwise.
     *                    If transa == 'N' or 'n', lda must be at least max(1, m);
     *                    otherwise, lda must be at least max(1, k).
     * @param[in] lda     leading dimension of two-dimensional array used to store matrix A.
     * @param[in] B       array of dimensions (ldb, n)
     *                    if transb == 'N' or 'n', and of dimensions (ldb, k) otherwise.
     *                    If transb == 'N' or 'n', ldb must be at least max(1, k);
     *                    otherwise, lda must be at least max(1, n).
     * @param[in] ldb     leading dimension of two-dimensional array used to store matrix B.
     * @param[in] beta    scalar multiplier applied to C.
     *                    If zero, C does not have to be a valid input
     * @param[in] C       array of dimensions (ldc,n);
     *                    ldc must be at least max(1,m).
     * @param[in] ldc     leading dimension of two-dimensional array used to store matrix C.
     * TODO[doxy] Is the following description correct?
     * @param[out] syncToken updated based on C = alpha * op(A) * op(B) + beta * C
     */
    template<typename ValueType>
    static void gemm(
        const CBLAS_ORDER order,
        const CBLAS_TRANSPOSE transa,
        const CBLAS_TRANSPOSE transb,
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
        const IndexType ldc,
        tasking::SyncToken* syncToken );

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
     * TODO[doxy] Is the following description correct?
     * @param[out] syncToken  contains the solution matrix X satisfying
     *                              op(A) * X = alpha * B
     *                           or X * op(A) = alpha * B
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
        const IndexType ldb,
        tasking::SyncToken* syncToken );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    static    bool initialized; //!< static initialization used for registration

    static bool registerInterface();//!< registration
};

}
/* namespace lama */

#endif // LAMA_CUDABLAS3_HPP_
