/**
 * @file OpenMPLAPACK.hpp
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
 * @brief OpenMPLAPACK.hpp
 * @author lschubert
 * @date 02.07.2012
 * $Id$
 */
#ifndef LAMA_OPENMPLAPACK_HPP_
#define LAMA_OPENMPLAPACK_HPP_

// for dll_import
#include <lama/config.hpp>

// others
#include <lama/LAMATypes.hpp>
#include <lama/SyncToken.hpp>

#include <lama/openmp/BLASHelper.hpp>

// logging
#include <logging/logging.hpp>

#include <omp.h>

namespace lama
{

class OpenMPLAPACK
{
public:

    /**
     * @brief getrf_cpu computes the LU factorization of a general m-by-n
     * matrix A in floating point single precision as
     *      A = P*L*U,
     * where P is a permutation matrix, L is lower triangular with unit diagonal
     * elements (lower trapezoidal if m > n) and U is upper triangular (upper
     * trapezoidal if m < n). The routine uses partial pivoting, with row
     * interchanges. L and U will be stored within A whereby the diagonal elements
     * of L will not be stored.
     *
     * @param[in] order   Specifies, whether the matrix is stored in column major
     *                     order (i.e. CblasColMajor) or in row major order (i.e.
     *                     CblasRowMajor).
     * @param[in] m         Number of rows of matrix A; m must be at least zero. m
     *                       specifies, how many rows of A will be touched by this
     *                       function.
     * @param[in] n         Number of columns of matrix A; n must be at least zero.
     *                       n specifies, how many columns of A will be touched by
     *                       this function.
     * @param[in,out] a     On input, a is holding the array, containing the matrix
     *                       A. On output it will be overwritten by L and U, whereby
     *                       the diagonal elements of L will not be stored.
     * @param[in] lda       The first dimension of array a. lda specifies the actual
     *                       number of rows of A. lda is used to compute the
     *                       position of the next column, i.e. if position (r,c) is
     *                       going to be accessed its position will be computed by
     *                       (c*lda + r).
     * @param[in,out] ipiv  The array holding the permutation of the matrix. Its
     *                       size must be at least min(m,n). It contains the info
     *                       that row i has been changed with ipiv[i].
     *                      ipiv is assumed to be stored in C-Style, i.e the values,
     *                       representing the indexes of the matrix are assumed to
     *                       be starting with zero and ending with m-1. It will
     *                       also leave with this assumption, but in between, it is
     *                       first incremented and afterwards decremented, to fit
     *                       the Fortran interface.
     *
     * @return info         If info=0, the execution is successful.
     *                      If info = -i, the i-th parameter had an illegal value.
     *                      If info = i, uii is 0. The factorization has been
     *                       completed, but U is exactly singular. Division by 0
     *                       will occur if you use the factor U for solving a
     *                       system of linear equations.
     */
    template<typename T>
    static IndexType getrf(
        const enum CBLAS_ORDER order,
        const IndexType m,
        const IndexType n,
        T* const a,
        const IndexType lda,
        IndexType* const ipiv );

    /**
     *
     */
    template<typename T>
    static IndexType getri(
        const enum CBLAS_ORDER order,
        const IndexType n,
        T* const A,
        const IndexType lda,
        IndexType* const ipiv );

    /** Implementation of LAPACKInterface::getinv vi LAPACK. */

    template<typename T>
    static void getinv( const IndexType n, T* a, const IndexType lda );

    /**
     * @brief trtrs solves the following equation system:
     *      op(A)*X = B
     *  where op(A) is either A, AT or AH;
     *  and B is a matrix of right hand sides and will contain the solution of all
     *  equations on output.
     *
     * @param[in] order   Specifies, whether the matrix is stored in column major
     *                     order (i.e. CblasColMajor) or in row major order (i.e.
     *                     CblasRowMajor).
     * @param[in] uplo    Specifies, whether matrix A is upper triangular (i.e.
     *                     CblasUpper) or lower triangular (i.e. CblasLower).
     * @param[in] trans   Specifies op(A).
     *                     if trans == CblasNoTrans,   op(A) = A;
     *                     if trans == CblasTrans,     op(A) = AT;
     *                     if trans == CblasConjTrans, op(A) = AH;
     * @param[in] diag    Specifies, whether the triangualr matrix A is a unit
     *                     triangular matrix, i.e. the diagonal elements of A are
     *                     one.
     *                     if diag == CblasNonUnit, the diagonal elements of A are
     *                      not assumed to be one;
     *                     if diag == CblasUnit, the diagonal elements of A are
     *                      assumed to be one and therefor not referenced;
     * @param[in] n       The number of columns of A and rows of B; n must be at
     *                     least 0.
     * @param[in] nrhs    The number of columns of B; nrhs must be at least 0.
     * @param[in] A       The array containing matrix A.
     * @param[in] lda     The first dimension of A, i.e. the actual number of rows
     *                     of A.
     * @param[in,out] B   On input B is the array holding the right hand sides of
     *                     the equations, on output, it will hold the solution for
     *                     each equation system.
     * @param[in] ldb     The first dimension of B, i.e. the actual number of columns
     *                     of B.
     */
    template<typename T>
    static IndexType trtrs(
        const enum CBLAS_ORDER order,
        const enum CBLAS_UPLO uplo,
        const enum CBLAS_TRANSPOSE trans,
        const enum CBLAS_DIAG diag,
        const IndexType n,
        const IndexType nrhs,
        const T* A,
        const IndexType lda,
        T* B,
        const IndexType ldb );

    /**
     * @brief tptrs solves the following equation system:
     *      op(A)*X = B
     *  where op(A) is either A, AT or AH;
     *  and B is a matrix of right hand sides and will contain the solution of all
     *  equations on output.
     *
     * @param[in] order   Specifies, whether the matrix is stored in column major
     *                     order (i.e. CblasColMajor) or in row major order (i.e.
     *                     CblasRowMajor).
     * @param[in] uplo    Specifies, whether matrix A is upper triangular (i.e.
     *                     CblasUpper) or lower triangular (i.e. CblasLower).
     * @param[in] trans   Specifies op(A).
     *                     if trans == CblasNoTrans,   op(A) = A;
     *                     if trans == CblasTrans,     op(A) = AT;
     *                     if trans == CblasConjTrans, op(A) = AH;
     * @param[in] diag    Specifies, whether the triangualr matrix A is a unit
     *                     triangular matrix, i.e. the diagonal elements of A are
     *                     one.
     *                     if diag == CblasNonUnit, the diagonal elements of A are
     *                      not assumed to be one;
     *                     if diag == CblasUnit, the diagonal elements of A are
     *                      assumed to be one and therefor not referenced;
     * @param[in] n       The number of columns of A and rows of B; n must be at
     *                     least 0.
     * @param[in] nrhs    The number of columns of B; nrhs must be at least 0.
     * @param[in] AP      The array containing matrix A.
     * @param[in,out] B   On input B is the array holding the right hand sides of
     *                     the equations, on output, it will hold the solution for
     *                     each equation system.
     * @param[in] ldb     The first dimension of B, i.e. the actual number of columns
     *                     of B.
     */
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

    /**
     * @brief laswp_cpu performs a series of row interchanges on the matrix A.
     * One row interchange is initiated for each of rows k1 through k2 of A.
     *
     * @param[in] order      Specifies, whether the matrix is stored in column major
     *                       order (i.e. CblasColMajor) or in row major order (i.e.
     *                       CblasRowMajor). Since a translation of the data would be
     *                       too expensiv, if it was stored in row major order, the
     *                       BLAS level1 function SSWAP will be called instead. The
     *                       beginning column of the vector in A will then be LDA-N.
     * @param[in] n          The number of columns of the matrix A.
     * @param[in,out] A      Array of dimension (LDA,N). On entry, the matrix of
     *                       column dimension N to which the row interchanges will be
     *                       applied. On exit, the permuted matrix.
     * @param[in] lda        If the matrix is stored in column major order, lda
     *                       specifies the actual number of rows of A. If else the
     *                       matrix is stored in row major order, lda specifies the
     *                       actual number of columns of A.
     * @param[in] k1         The first element of ipiv for which a row interchange will
     *                       be done.
     * @param[in] k2         The last element of ipiv for which a row interchange will
     *                       be done.
     * @param[in] ipiv       Array of dimension (k2*abs(incx)). The vector of pivot
     *                       indices. Only the elements in positions k1 through k2 of
     *                       ipiv are accessed. ipiv(k) = l implies rows k and l are
     *                       to be interchanged.
     * @param[in] incx       The increment between successive values of ipiv. If ipiv
     *                       is negative, the pivots are applied in reverse order.
     * @param[out] syncToken TODO[doxy] Complete Description.
     */
    template<typename T>
    static void laswp(
        const enum CBLAS_ORDER order,
        const IndexType n,
        T* A,
        const IndexType lda,
        const IndexType k1,
        const IndexType k2,
        const IndexType* ipiv,
        const IndexType incx,
        SyncToken* syncToken );

    /** Routine that sets functions pointers belonging to LAPACK in a BLASInterface.
     *
     *  param[inout] BLASInterface struct to register all routines implemented in OpenMP
     *
     *  Note: this routine will make instantiations of the template routines.
     */

    static void setInterface( struct BLASInterface& BLAS );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

}; /* OpenMPLAPACK */

} /* namespace lama */

#endif // LAMA_OPENMPLAPACK_HPP_
