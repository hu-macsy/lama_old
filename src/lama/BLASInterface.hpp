/**
 * @file BLASInterface.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Interface class for context dependent BLAS operations used in LAMA.
 * @author Thomas Brandes
 * @date 02.04.2013
 * $Id$
 */
#ifndef LAMA_BLAS_INTERFACE_HPP_
#define LAMA_BLAS_INTERFACE_HPP_

#include <lama/LAMATypes.hpp>

// others
#include <lama/macros/interface.hpp>

// C++
#include <cstring>

namespace lama
{

/** Interface struct for BLAS routines.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 *
 *  The interface is grouped in several structures 
 *
 *  - BLAS1 for BLAS level-1 routines
 *  - BLAS2 for BLAS level-2 routines
 *  - BLAS3 for BLAS level-3 routines
 *  - LAPACK for LAPACK routines
 *  - SCALAPACK for SCALAPACK routines
 *
 */

struct BLASInterface
{

    /** Function pointer type definitions for BLAS1 routines. */

    template<typename ValueType> 
    struct BLAS1
    {
        /**
         * @brief scal replaces vector x with alpha * x.
         *
         *        x = alpha * x
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] alpha  scalar multiplier
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * @param[out] x     vector x(unchanged if n<=0 or incX <=0)
         */
        typedef void ( *scal ) ( const IndexType n, 
                                 const ValueType alpha, 
                                 ValueType* x, 
                                 const IndexType incX, 
                                 class SyncToken* syncToken );

        /**
         * @brief nrm2 computes the Euclidean norm of the n-vector x
         * (with storage increment incX).
         *
         *         nrm2(x) = sqrt( sum_{i = 1}^{n}( x_i^2 ) )
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * return            the Euclidian norm (returns zero if n <=0, incX <=0)
         */
        typedef ValueType ( *nrm2 ) ( const IndexType n, 
                                      const ValueType* x, 
                                      const IndexType incX, 
                                      class SyncToken* syncToken );

        /**
         * @brief asum computes the sum of the absolute values of the elements
         * of vector x
         *
         *          asum(x) = sum_{i = 1}^{n}( x_i )
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * return            the sum of absolute values (returns zero if n<=0 or incX <=0)
         */
        typedef ValueType ( *asum ) ( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken );
    
        /** max norm: maxnorm(x) = max( x_i ) */
    
        /**
         * @brief iamax finds the smallest index of the maximum magnitude
         * element of vector x
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * return            the smallest index (returns zero if n <=0 or incX <=0)
         */
        typedef IndexType (*iamax) ( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken );
    
        /**
         * @brief iamax finds the smallest index of the maximum magnitude
         * element of vector x
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * return            the smallest index (returns zero if n <=0 or incX <=0)
         */
        typedef ValueType (*viamax) ( const IndexType n, const ValueType* x, const IndexType incX, SyncToken* syncToken );
    
        /**
         * @brief swap interchanges vector x with vector y.
         *
         *            x <-> y
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * @param[in] y      vector with n elements
         * @param[in] incY   storage spacing between elements of y
         * @param[out] x     vector x(unchanged if n<=0)
         * @param[out] y     vector y(unchanged if n<=0)
         */
        typedef void (*swap) ( const IndexType n, ValueType* x, const IndexType incX, ValueType* y, const IndexType incY, SyncToken* syncToken );
    
        /**
         * @brief copy copies the vector x to the vector y.
         *
         *            y = x
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * @param[in] y      vector with n elements
         * @param[in] incY   storage spacing between elements of y
         * @param[out] y     contains vector x
         */
        typedef void ( *copy )( const IndexType n,
                                const ValueType* x,
                                const IndexType incX,
                                ValueType* y,
                                const IndexType incY,
                                class SyncToken* syncToken );
    
        /**
         * @brief axpy multiplies vector x by scalar alpha and
         * adds the result to vector y.
         *
         *            y = alpha * y + beta * x
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] alpha  scalar multiplier
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * @param[in] y      vector with n elements
         * @param[in] incY    storage spacing between elements of y
         * @param[out] y     result (unchanged if n<=0)
         */
        typedef void ( *axpy )( const IndexType n,
                                ValueType alpha,
                                const ValueType* x,
                                const IndexType incX,
                                ValueType* y,
                                const IndexType incY,
                                class SyncToken* syncToken );
    
        /**
         * @brief dot computes the dot product of two vectors.
         * It returns the dot product of the vectors x and y if successfull,
         * and 0.0 otherwise.
         * adds the result to vector y.
         *
         *            dot = sum_{i = 1}^{n}( x_i * y_i )
         *
         * @param[in] n      number of elements in input vectors.
         * @param[in] x      vector with n elements
         * @param[in] incX   storage spacing between elements of x
         * @param[in] y      vector with n elements
         * @param[in] incY   storage spacing between elements of y
         * return            dot product (returns zero if n <= 0)
         */
        typedef ValueType ( *dot )( const IndexType n,
                                    const ValueType* x,
                                    const IndexType incX,
                                    const ValueType* y,
                                    const IndexType inc,
                                    class SyncToken* syncToken );
    
        /**  sum: z = alpha * x + beta * y */
        typedef void ( *sum ) ( const IndexType n, 
                                const ValueType alpha, 
                                const ValueType* x, 
                                ValueType beta, 
                                const ValueType* y, 
                                ValueType* z, 
                                class SyncToken* syncToken );
    
        /**
         * @brief ass The function ass() assigns one scalar value to a vector of the given size.
         *
         * @param[in] n         size of the vector
         * @param[in] value     scalar value, which should be assign to the whole vector
         * @param[out] x        vector, the values should be assigned to
         */
         typedef void ( *ass ) ( const IndexType n, const ValueType value, ValueType *x, class SyncToken* syncToken );
    
    };

    // declare variables of function pointers, i.e. arrays indexed by each type

    LAMA_INTERFACE_DEFINE_T( BLAS1, scal )
    LAMA_INTERFACE_DEFINE_T( BLAS1, nrm2 )
    LAMA_INTERFACE_DEFINE_T( BLAS1, asum )
    LAMA_INTERFACE_DEFINE_T( BLAS1, iamax )
    LAMA_INTERFACE_DEFINE_T( BLAS1, viamax )
    LAMA_INTERFACE_DEFINE_T( BLAS1, copy )
    LAMA_INTERFACE_DEFINE_T( BLAS1, axpy )
    LAMA_INTERFACE_DEFINE_T( BLAS1, dot )
    LAMA_INTERFACE_DEFINE_T( BLAS1, sum )
    LAMA_INTERFACE_DEFINE_T( BLAS1, ass )

    template<typename ValueType>
    struct BLAS2
    {
        /**
         * @brief gemv performs one of the matrix-vector operations
         *
         * y = alpha * op(A) * x + beta * y
         * where op(A) = A or op(A) = AT
         *
         * alpha and beta are scalars, and x and y are vectors.
         * A is an m×n matrix consisting of elements.
         * Matrix A is stored in column-major format, and lda is the
         * leading dimension of the two dimensional array in which A is stored.
         *
         * @param[in] trans    specifies op(A). If trans == 'N' or 'n', op(A) = A
         *                                      If trans == 'T','t','C','c', op(A) = AT
         * @param[in] m        the number of rows of matrix A; m must be at least zero
         * @param[in] n        the number of columns of matrix A; n must be at least zero
         * @param[in] alpha    scalar multiplier applied to op(A)
         * @param[in] A        array of dimensions (lda,n).
         *                     If trans == 'N' or 'n', of dimensions (lda,m)
         *                     otherwise; lda must be at least max(1,m) if trans == 'N' or 'n'
         *                     and at least max(1,n) otherwise
         * @param[in] lda      leading dimension of two-dimensional array used to store matrix A.
         * @param[in] x        array of length at least (1 + (n - 1) * abs(incX))
         *                     when trans == 'N' or 'n'
         *                     else at least (1 + ( m - 1) * abs(incX))
         * @param[in] incX      storage spacing between elements of x; incX must not be zero
         * @param[in] beta      scalar multiplier applied to vector y. If beta is zero, y is not read
         * @param[in] y         array of length at least (1 + (m - 1) * abs(incY))
         *                      when trans == 'N' or 'n'
         *                      else at least (1 + ( n - 1) * abs(incY))
         * @param[in] incY      the storage spacing between elements of y; incY must not be zero.
         * @param[out] y        updated according to y = alpha * op(A) * x + beta * y
         * @param[in] syncToken allows to start asynchronous execution
         *
         */
        typedef void ( *gemv ) ( const enum CBLAS_ORDER order,
                                 const enum CBLAS_TRANSPOSE trans,
                                 const IndexType m,
                                 const IndexType n,
                                 const ValueType alpha,
                                 const ValueType *A,
                                 const IndexType lda,
                                 const ValueType *x,
                                 const IndexType incX,
                                 const ValueType beta,
                                 ValueType *y,
                                 const IndexType incY,
                                 class SyncToken* syncToken );
    };
    
    LAMA_INTERFACE_DEFINE_T( BLAS2, gemv )

    template<typename ValueType>
    struct BLAS3
    {
        /**
         * @brief gemm computes the product of matrix A and matrix B,
         * multiplies the result by scalar alpha,
         * and adds the sum to the product of matrix C and scalar beta.
         * It performs one of the matrix‐matrix operations:
         *
         * C = alpha * op(A) * op(B) + beta * C,
         * where op(X) = X or op(X) = \f$X^{T}\f$
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
         * @param[in] order   specifies whether arrays are row-major or column-major stored
         * @param[in] transA  specifies op(A)
         *                    If transa = 'N', 'n'
         *                    op(A) = A
         *                    If transa == 'T', 't', 'C', 'c'
         *                    op(A) = \f$A^{T}\f$
         * @param[in] transB  specifies op(B)
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
         * @param[in,out] C   array of dimensions (ldc,n); ldc must be at least max(1,m).
         *                    updated based on C = alpha * op(A) * op(B) + beta * C
         * @param[in] ldc     leading dimension of two-dimensional array used to store matrix C.
         * @param[in,out] syncToken is optional synchronization taken that might be used for asynchronous execution
         */
    
        typedef void ( *gemm ) ( const enum CBLAS_ORDER order, 
                                 const enum CBLAS_TRANSPOSE transA, 
                                 const enum CBLAS_TRANSPOSE transB,
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
                                 class SyncToken* syncToken );
    
    };

    LAMA_INTERFACE_DEFINE_T( BLAS3, gemm )

    /** Structure with pointers for routines using LAPACK. */

    template<typename ValueType>
    struct LAPACK
    {
        /**  */
        typedef IndexType (*getrf) (const enum CBLAS_ORDER order, const IndexType m, const IndexType n, ValueType* a,
                            const IndexType lda, IndexType* ipivot);

        /** Method computes the inverse of a matrix by using the LAPACK routines getrf and getri
         *
         *  @param[in]     n specifies the order of the matrix a
         *  @param[in,out] a is the matrix for which the inverse is computed in-place
         *  @param[in]     lda for the leading dimension of the array A
         *  @throws        Exception if error occurs ( e.g. matrix is singular )
         *
         *  Note that the storage order (column-wise or row-wise does not matter at all)
         */
    
        typedef void (*getinv) ( const IndexType n, ValueType* a, const IndexType lda );
    
        /**  */
        typedef IndexType (*getri) (const enum CBLAS_ORDER , const IndexType n, ValueType* a,
                            const IndexType lda, IndexType* ipivot);

        /**  */
        typedef IndexType (*tptrs) (const enum CBLAS_ORDER order, const enum CBLAS_UPLO uplo, const enum CBLAS_TRANSPOSE trans,
                            const enum CBLAS_DIAG diag, const IndexType n, const IndexType nrhs, const ValueType* AP, ValueType* B,
                            const IndexType ldb);
    
        /**  */
        typedef void ( *laswp ) ( const enum CBLAS_ORDER order, 
                                  const IndexType n, 
                                  ValueType* A, 
                                  const IndexType lda, 
                                  const IndexType k1,
                                  const IndexType k2,   
                                  const IndexType* ipiv,   
                                  const IndexType incx,   
                                  class SyncToken* syncToken);
    
    };

    LAMA_INTERFACE_DEFINE_T( LAPACK, getrf )
    LAMA_INTERFACE_DEFINE_T( LAPACK, getri )
    LAMA_INTERFACE_DEFINE_T( LAPACK, getinv )
    LAMA_INTERFACE_DEFINE_T( LAPACK, trtrs )
    LAMA_INTERFACE_DEFINE_T( LAPACK, tptrs )
    LAMA_INTERFACE_DEFINE_T( LAPACK, laswp )

    template<typename ValueType>
    struct SCALAPACK
    {
        /** Function pointer for routine that computes the inverse of a cyclic(nB) distributed matrix.
         *
         *  @param[in]  n  global size of the matrix,
         *  @param[in]  a  is pointer to the values of the local dense storage
         *  @param[in]  nb is the blocking factor of the cyclic distribution
         *  @param[in]  comm is the communicator of the distribution
         */

        typedef void ( *inverse ) ( const IndexType n, const IndexType nB, const ValueType* a, const class Communicator& comm );
    };
    
    LAMA_INTERFACE_DEFINE_T( SCALAPACK, inverse )

    /** Default constructor, initializes function pointer variables with NULL */

    BLASInterface()
    {
        // intialize all function pointers with NULL

        memset( this, 0, sizeof( *this ) );
    }
};

} //namespace lama

#endif // LAMA_BLAS_INTERFACE_HPP_
