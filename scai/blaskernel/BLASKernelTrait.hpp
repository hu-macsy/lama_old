/**
 * @file BLASKernelTrait.hpp
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
 * @brief Kernel traits for BLAS and LAPACK operations used in LAMA.
 * @author Thomas Brandes
 * @date 02.04.2013
 */
#pragma once

// local libary
#include <scai/blaskernel/cblas.hpp>

// internal scai libraries
#include <scai/common/SCAITypes.hpp>
#include <scai/common/MatrixOp.hpp>

// std
#include <cstring>

namespace scai
{

/** Namespace for all kernel implementations of BLAS routines. */

namespace blaskernel
{

/** Struct with all LAMA kernel traits for BLAS and LAPACK routines currently used.
 *
 *  Syntax of the LAMA kernel definitions are nearly the same as used for the
 *  corresponing BLAS and LAPACK routines.
 *
 *  It could have been just a namespace but struct disables
 *  'using namespace BLASKerneltrait'.
 */

struct BLASKernelTrait
{

    /** Kernel trait for BLAS1 routine scal.
     *
     *  @tparam ValueType stands for the arithmetic type used for scal
     */
    template<typename ValueType>
    struct scal
    {
        /**
         * @brief scal replaces vector x with alpha * x.
         *
         *        x = alpha * x
         *
         * @param[in] n      number of considered elements in input vector.
         * @param[in] alpha  scalar multiplier.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.
         * @param[in] incX   incX storage spacing between considered elements of x.
         * @param[out] x     vector x (unchanged if n<=0 or incX <= 0 ).
         */
        typedef void ( *FuncType )(
            const IndexType n,
            const ValueType alpha,
            ValueType* x,
            const IndexType incX );

        static const char* getId()
        {
            return "BLAS1.scal";
        }
    };

    /** Kernel trait for BLAS routine nrm2.
     *
     *  @tparam ValueType stands for the arithmetic type used for nrm2
     */
    template<typename ValueType>
    struct nrm2
    {
        /**
         * @brief nrm2 computes the Euclidean norm of the n-vector x
         * (with storage increment incX).
         *
         *         nrm2(x) = sqrt( sum_{i = 1}^{n}( x_i^2 ) )
         *
         * @param[in] n      number of considered elements in input vector.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.
         * @param[in] incX   incX storage spacing between considered elements of x.
         * return            the Euclidian norm within considered elements
         *                   of x (returns zero if n <=0, incX <= 0).
         */

        typedef ValueType ( *FuncType )( const IndexType n, const ValueType* x, const IndexType incX );

        static const char* getId()
        {
            return "BLAS1.nrm2";
        }
    };

    /** Kernel trait for BLAS routine asum.
     *
     *  @tparam ValueType stands for the arithmetic type used for asum
     */
    template<typename ValueType>
    struct asum
    {

        /**
         * @brief asum computes the sum of the absolute values of the elements
         * of vector x
         *
         *          asum(x) = sum_{i = 1}^{n}( x_i )
         *
         * @param[in] n      number of considered elements in input vectors.
         * @param[in] x      vector with minimum (n - 1) * incX + 1  elements.
         * @param[in] incX   incX storage spacing between considered elements of x.
         * return            the sum of absolute values within considered elements
         *                   of x (returns zero if n<=0 or incX <= 0).
         */
        typedef ValueType ( *FuncType )( const IndexType n, const ValueType* x, const IndexType incX );

        static const char* getId()
        {
            return "BLAS1.asum";
        }
    };

    /** Kernel trait for BLAS routine iamax.
     *
     *  @tparam ValueType stands for the arithmetic type used for iamax
     */
    template<typename ValueType>
    struct iamax
    {

        /** max norm: maxnorm(x) = max( x_i ) */

        /**
         * @brief iamax finds the smallest index of the maximum magnitude
         * element of vector x.
         *
         * @param[in] n      number of considered elements in input vector.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.iama
         *
         * @param[in] incX   storage spacing between considered elements of x.
         * return            the smallest index of the maximum magnitude
         *                   within considered elements of x (returns zero if n <=0 or incX <=0).
         */
        typedef IndexType ( *FuncType )( const IndexType n, const ValueType* x, const IndexType incX );

        static const char* getId()
        {
            return "BLAS1.iamax";
        }
    };

    /** Kernel trait for BLAS1 routine swap.
     *
     *  @tparam ValueType stands for the arithmetic type used for iamax
     */
    template<typename ValueType>
    struct swap
    {

        /**
         * @brief swap interchanges considered elements of vector x with vector y.
         *
         *            x <-> y
         *
         * @param[in] n      number of considered elements in input vectors.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.
         * @param[in] incX   storage spacing between considered elements of x.
         * @param[in] y      vector with minimum (n - 1) * incY + 1 elements.
         * @param[in] incY   storage spacing between considered elements of y.
         * @param[out] x     vector x (unchanged if n<=0, incX <=0 or or incY <=0).
         * @param[out] y     vector y (unchanged if n<=0, incX <=0 or or incY <=0).
         */
        typedef void ( *FuncType )(
            const IndexType n,
            ValueType* x,
            const IndexType incX,
            ValueType* y,
            const IndexType incY );

        static const char* getId()
        {
            return "BLAS1.copy";
        }
    };

    /** Kernel trait for BLAS1 routine copy.
     *
     *  @tparam ValueType stands for the arithmetic type used for copy
     */
    template<typename ValueType>
    struct copy
    {

        /**
         * @brief copy copies the vector x to the vector y.
         *
         *            y = x
         *
         * @param[in] n      number of considered elements in input vectors.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements
         * @param[in] incX   storage spacing between elements of x
         * @param[in] y      vector with minimum (n - 1) * incY + 1 elements
         * @param[in] incY   storage spacing between elements of y
         * @param[out] y     result (unchanged if n<=0, incX <=0 or or incY <=0 )
         */
        typedef void ( *FuncType )(
            const IndexType n,
            const ValueType* x,
            const IndexType incX,
            ValueType* y,
            const IndexType incY );

        static const char* getId()
        {
            return "BLAS1.copy";
        }
    };

    /** Kernel trait for BLAS1 routine axpy.
     *
     *  @tparam ValueType stands for the arithmetic type used for axpy.
     */
    template<typename ValueType>
    struct axpy
    {

        /**
         * @brief axpy multiplies scalar alpha by vector x and
         * adds the result to vector y.
         *
         *            y = y + alpha * x
         *
         * @param[in] n      number of considered elements in input vectors.
         * @param[in] alpha  scalar multiplier
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.
         * @param[in] incX   storage spacing between elements of x
         * @param[in,out] y  vector with minimum (n - 1) * incY + 1 elements.
         * @param[in] incY   storage spacing between elements of y
         */
        typedef void ( *FuncType )(
            const IndexType n,
            ValueType alpha,
            const ValueType* x,
            const IndexType incX,
            ValueType* y,
            const IndexType incY );

        static const char* getId()
        {
            return "BLAS1.axpy";
        }
    };

    /** Kernel trait for BLAS1 routine dot.
     *
     *  @tparam ValueType stands for the arithmetic type used for dot.
     */
    template<typename ValueType>
    struct dot
    {

        /**
         * @brief dot computes the dot product of two vectors.
         * It returns the dot product of the vectors x and y if successfull,
         * and 0.0 otherwise.
         *
         *            dot = sum_{i = 0}^{n-1}( x_i * y_i )
         *
         * If ValueType is complex, the conj( x_i ) is taken
         *
         * @param[in] n      number of considered elements in input vectors.
         * @param[in] x      vector with minimum (n - 1) * incX + 1 elements.
         * @param[in] incX   storage spacing between elements of x.
         * @param[in] y      vector with minimum (n - 1) * incY + 1 elements.
         * @param[in] incY   storage spacing between elements of y.
         * @return           dot product (returns zero if n <= 0, incX <=0 or or incY <=0).
         */
        typedef ValueType ( *FuncType )(
            const IndexType n,
            const ValueType* x,
            const IndexType incX,
            const ValueType* y,
            const IndexType inc );

        static const char* getId()
        {
            return "BLAS1.dot";
        }
    };

    /** Kernel trait for LAMA specific sum routine.
     *
     *  @tparam ValueType stands for the arithmetic type used for scal
     */
    template<typename ValueType>
    struct sum
    {

        /**
         * @brief sum adds the multiplication of scalar alpha by vector x and the multiplication of
         * scalar beta by vector y into a new vector z.
         *
         *            z = alpha * x + beta * y
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] alpha  scalar multiplier.
         * @param[in] x      vector with n elements.
         * @param[in] beta   scalar multiplier.
         * @param[in] y      vector with n elements.
         * @param[out] z     result of adding the 2 multiplications (unchanged if n <= 0).
         */
        typedef void ( *FuncType )(
            const IndexType n,
            const ValueType alpha,
            const ValueType* x,
            ValueType beta,
            const ValueType* y,
            ValueType* z );

        static const char* getId()
        {
            return "BLAS1.sum";
        }
    };

    /** Kernel trait for BLAS2 routine gemv.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct gemv
    {
        /**
         *
         * y = alpha * op(A) * x + beta * y
         * where op(A) = A or op(A) = AT
         *
         * alpha and beta are scalars, and x and y are vectors.
         * A is an m×n matrix consisting of elements.
         * Matrix A is stored in row-major format, and lda is the
         * leading dimension of the two dimensional array in which A is stored.
         *
         * @param[in] opA       specifies op(A).
         * @param[in] m         the number of rows of matrix A; m must be at least zero
         * @param[in] n         the number of columns of matrix A; n must be at least zero
         * @param[in] alpha     scalar multiplier applied to op(A)
         * @param[in] A         array of dimensions (lda,n).
         *                      If trans == 'N' or 'n', of dimensions (lda,m)
         *                      otherwise; lda must be at least max(1,m) if trans == 'N' or 'n'
         *                      and at least max(1,n) otherwise
         * @param[in] lda       leading dimension of two-dimensional array used to store matrix A.
         * @param[in] x         array of length at least (1 + (n - 1) * incX)
         *                      when trans == 'N' or 'n'
         *                      else at least (1 + ( m - 1) * incX)
         * @param[in] incX      storage spacing between elements of x; incX must be >= 0.
         * @param[in] beta      scalar multiplier applied to vector y. If beta is zero, y is not read
         * @param[in] y         array of length at least (1 + (m - 1) * incY)
         *                      when trans == 'N' or 'n'
         *                      else at least (1 + ( n - 1) * incY)
         * @param[in] incY      the storage spacing between elements of y; incY must be >= 0.
         * @param[out] y        updated according to y = alpha * op(A) * x + beta * y
         *
         */
        typedef void ( *FuncType ) (
            const common::MatrixOp opA,
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

        static const char* getId()
        {
            return "BLAS2.gemv";
        }
    };

    /** Kernel trait for BLAS2 extension routine matrix add
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct geam
    {
        typedef void ( *FuncType ) (
            ValueType* C,
            const IndexType ldc,
            const IndexType m,
            const IndexType n,
            const ValueType alpha,
            const ValueType* A,
            const IndexType lda,
            const common::MatrixOp opA,
            const ValueType beta,
            const ValueType* B,
            const IndexType ldb,
            const common::MatrixOp opB );

        static const char* getId()
        {
            return "BLAS2.geam";
        }
    };

    /** Kernel trait for BLAS3 routine gemm.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct gemm
    {
        /**
         * @brief gemm computes the product of matrix A and matrix B,
         * multiplies the result by scalar alpha,
         * and adds the sum to the product of matrix C and scalar beta.
         * It performs one of the matrix‐matrix operations:
         *
         * C = alpha * op(A) * op(B) + beta * C,
         * where op is one of the supported matrix operations conj/transpose
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
         * @param[in] opA     specifies op(A), e.g. NORMAL, TRANSPOSE, CONJ, CONJ_TRANSPOSE
         * @param[in] opB     specifies op(B), e.g. NORMAL, TRANSPOSE, CONJ, CONJ_TRANSPOSE
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
         * @param[in,out] C   array of dimensions (ldc,n); ldc must be at least max(1,m).
         *                    updated based on C = alpha * op(A) * op(B) + beta * C
         * @param[in] ldc     leading dimension of two-dimensional array used to store matrix C.
         *
         * All matrices are assumed to be in the row-major format.
         */

        typedef void ( *FuncType ) (
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

        static const char* getId()
        {
            return "BLAS3.gemm";
        }
    };

    /** Kernel trait for LAPACK routine getrf.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct getrf
    {
        /* @brief computes the LU factorization of a general m-by-n
         * matrix A in floating point single precision as
         *      A = P*L*U,
         * where P is a permutation matrix, L is lower triangular with unit diagonal
         * elements (lower trapezoidal if m > n) and U is upper triangular (upper
         * trapezoidal if m < n). The routine uses partial pivoting, with row
         * interchanges. L and U will be stored within A whereby the diagonal elements
         * of L will not be stored.
         *
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
         * @throws common::Exception for illegal arguments or incomplete factorization
         */

        typedef void ( *FuncType ) (
            const IndexType m,
            const IndexType n,
            ValueType* a,
            const IndexType lda,
            IndexType* ipivot );

        static const char* getId()
        {
            return "LAPACK.getrf";
        }
    };

    /** Kernel trait for getinv as concatenation of getrf and getri.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct getinv
    {
        /** Method computes the inverse of a matrix by using the LAPACK routines getrf and getri
         *
         *  @param[in]     n specifies the size of the matrix a in each dimension
         *  @param[in,out] a is the matrix for which the inverse is computed in-place
         *  @param[in]     lda for the leading dimension of the array A
         *  @throws        Exception if error occurs ( e.g. matrix is singular )
         *
         *  Note that the storage order (column-wise or row-wise does not matter at all)
         */

        typedef void ( *FuncType ) ( const IndexType n, ValueType* a, const IndexType lda );

        static const char* getId()
        {
            return "LAPACK.getinv";
        }
    };

    /** Kernel trait for LAPACK routine getri.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct getri
    {
        typedef void ( *FuncType ) (
            const IndexType n,
            ValueType* a,
            const IndexType lda,
            IndexType* ipivot );

        static const char* getId()
        {
            return "LAPACK.getri";
        }
    };

    /** Kernel trait for LAPACK routine tptrs.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct tptrs
    {
        /**
         * @brief tptrs solves the following equation system:
         *      op(A)*X = B
         *  where op(A) is either A, AT or AH;
         *  and B is a matrix of right hand sides and will contain the solution of all
         *  equations on output.

         * @param[in] uplo    Specifies, whether matrix A is upper triangular (i.e.
         *                    CblasUpper) or lower triangular (i.e. CblasLower).
         * @param[in] op      Specifies op(A).
         *                    if op == MatrixOp::NORMAL,         op(A) = A;
         *                    if op == MatrixOp::TRANSPOSE,      op(A) = AT;
         *                    if op == MatrixOp::CONJ_TRANSPOSE, op(A) = AH;
         * @param[in] diag    Specifies, whether the triangualr matrix A is a unit
         *                    triangular matrix, i.e. the diagonal elements of A are
         *                    one.
         *                    if diag == CblasNonUnit, the diagonal elements of A are
         *                    not assumed to be one;
         *                    if diag == CblasUnit, the diagonal elements of A are
         *                    assumed to be one and therefor not referenced;
         * @param[in] n       number of columns of A and rows of B; n must be at least 0.
         * @param[in] nrhs    The number of columns of B; nrhs must be at least 0.
         * @param[in] AP      The array containing matrix A.
         * @param[in,out] B   On input B is the array holding the right hand sides of
         *                     the equations, on output, it will hold the solution for
         *                     each equation system.
         * @param[in] ldb     The first dimension of B, i.e. the actual number of columns
         *                     of B.
         */

        typedef void ( *FuncType ) (
            const CBLAS_UPLO uplo,
            const common::MatrixOp op,
            const CBLAS_DIAG diag,
            const IndexType n,
            const IndexType nrhs,
            const ValueType* AP,
            ValueType* B,
            const IndexType ldb );

        static const char* getId()
        {
            return "LAPACK.tptrs";
        }
    };

    /** Kernel trait for using SCALAPACK to compute inverse.
     *
     *  @tparam ValueType stands for the arithmetic type used in this operation.
     */
    template<typename ValueType>
    struct inverse
    {
        /** Function pointer for routine that computes the inverse of a cyclic(nB) distributed matrix.
         *  (here SCALAPACK where kernel itself might call MPI, so it must be enabled before)
         *
         *  @param[in]  n  global size of the matrix,
         *  @param[in]  a  is pointer to the values of the local dense storage
         *  @param[in]  nb is the blocking factor of the cyclic distribution
         *  @param[in]  comm is the communicator of the distribution
         */

        typedef void ( *FuncType ) ( const IndexType n, const IndexType nB, const ValueType* a, const class Communicator& comm );

        static const char* getId()
        {
            return "SCALAPACK.inverse";
        }
    };
};

} /* end namespace blaskernel */

} /* end namespace scai */
