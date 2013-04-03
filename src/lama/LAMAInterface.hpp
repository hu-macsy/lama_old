/**
 * @file LAMAInterface.hpp
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
 * @brief Interface class for context dependent operations to be implemented.
 * @author Thomas Brandes
 * @date 27.04.2011 (revised 02.04.2013)
 * $Id$
 */
#ifndef LAMA_LAMA_INTERFACE_HPP_
#define LAMA_LAMA_INTERFACE_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/Printable.hpp>

// others
#include <lama/Scalar.hpp>

/** Macro that creates for a function pointer type a table that contains pointer for each type.
 */

#define LAMA_INTERFACE_DEFINE( structname, functionname )                                     \
                                                                                              \
    /** Method that returns the function pointer for the correspondingly typed version. */    \
                                                                                              \
    structname::functionname                                                                  \
    functionname () const                                                                     \
    {                                                                                         \
        return functionname##_Table;                                                          \
    }                                                                                         \
                                                                                              \
    /** Register a function pointer in the type table. */                                     \
                                                                                              \
    void functionname##_add( structname::functionname functionPtr )                           \
    {                                                                                         \
        functionname##_Table = functionPtr;                                                   \
    }                                                                                         \
                                                                                              \
    structname::functionname functionname##_Table;

#define LAMA_INTERFACE_DEFINE_T( structname, functionname )                                   \
                                                                                              \
    /** Method that returns the function pointer for the correspondingly typed version. */    \
                                                                                              \
    template<typename T>                                                                      \
    typename structname<T>::functionname                                                      \
    functionname () const                                                                     \
    {                                                                                         \
        return ( typename structname<T>::functionname )                                       \
               functionname##_Table[ Scalar::getType<T>() ];                                  \
    }                                                                                         \
                                                                                              \
    /** Register a function pointer in the type table. */                                     \
                                                                                              \
    template<typename T>                                                                      \
    void functionname##_add( typename structname<T>::functionname functionPtr )               \
    {                                                                                         \
        functionname##_Table[ Scalar::getType<T>() ]                                          \
         = ( void (*) () ) functionPtr;                                                       \
    }                                                                                         \
                                                                                              \
    void ( *functionname##_Table[ Scalar::UNKNOWN ] ) ();

#define LAMA_INTERFACE_DEFINE_TT( structname, functionname )                                  \
                                                                                              \
    /** Method that returns the function pointer for the correspondingly typed version. */    \
                                                                                              \
    template<typename T1, typename T2>                                                        \
    typename structname<T1, T2>::functionname                                                 \
    functionname () const                                                                     \
    {                                                                                         \
        return ( typename structname<T1,T2>::functionname )                                   \
               functionname##_Table[ Scalar::getType<T1>()][ Scalar::getType<T2>() ];         \
    }                                                                                         \
                                                                                              \
    /** Register a function pointer in the type table. */                                     \
                                                                                              \
    template<typename T1, typename T2>                                                        \
    void functionname##_add( typename structname<T1, T2>::functionname functionPtr )          \
    {                                                                                         \
        functionname##_Table[ Scalar::getType<T1>() ][ Scalar::getType<T2>() ]                \
           = ( void (*) () ) functionPtr;                                                     \
    }                                                                                         \
                                                                                              \
    void ( *functionname##_Table[ Scalar::UNKNOWN ][ Scalar::UNKNOWN] ) ();

/** This macro intializes type table with of function will NULL function pointers.
 *  Macro can only be used within the interface itself
 */

#define LAMA_INTERFACE_INIT( functionname )                                                 \
    functionname##_Table = NULL;

#define LAMA_INTERFACE_INIT_T( functionname )                                               \
    for ( int i = 0; i < Scalar::UNKNOWN; ++i )                                             \
    {                                                                                       \
        functionname##_Table[i] = NULL;                                                     \
    }                                                                                       \

#define LAMA_INTERFACE_INIT_TT( functionname )                                              \
    for ( int i = 0; i < Scalar::UNKNOWN; ++i )                                             \
    {                                                                                       \
        for ( int j = 0; j < Scalar::UNKNOWN; ++j )                                         \
        {                                                                                   \
            functionname##_Table[i][j] = NULL;                                              \
        }                                                                                   \
    }                                                                                       \

#define LAMA_INTERFACE_REGISTER( interface, function )                                      \
    interface.function##_add( function );

#define LAMA_INTERFACE_REGISTER_T( interface, function, T )                                 \
    interface.function##_add<T>( function<T> );

#define LAMA_INTERFACE_REGISTER_TT( interface, function, T1, T2 )                           \
    interface.function##_add<T1,T2>( function<T1,T2> );

#define LAMA_INTERFACE_FN( function, loc, module, structname )                                \
    typename module##Interface::structname::function function =                               \
       loc->getInterface().module.function();                                                 \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

#define LAMA_INTERFACE_FN_DEFAULT( function, loc, module, structname )                        \
    typename module##Interface::structname::function function =                               \
       loc->getInterface().module.function();                                                 \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_UNSUPPORTED( "Method " #module "::" #function " not available on " << *loc );    \
        loc = ContextFactory::getContext( Context::Host );                                    \
        function = loc->getInterface().module.function();                                     \
        if ( function == NULL )                                                               \
        {                                                                                     \
            LAMA_THROWEXCEPTION( "Method " #module "::" #function                             \
                                 " also not available on " << *loc );                         \
        }                                                                                     \
    }

#define LAMA_INTERFACE_FN_T( function, loc, module, structname, T )                           \
    typename module##Interface::structname<T>::function function;                             \
    function = loc->getInterface().module.function<T>();                                      \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

#define LAMA_INTERFACE_FN_DEFAULT_T( function, loc, module, structname, T )                   \
    typename module##Interface::structname<T>::function function =                            \
       loc->getInterface().module.function<T>();                                              \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_UNSUPPORTED( "Method " #module "::" #function " not available on " << *loc );    \
        loc = ContextFactory::getContext( Context::Host );                                    \
        function = loc->getInterface().module.function<T>();                                  \
        if ( function == NULL )                                                               \
        {                                                                                     \
            LAMA_THROWEXCEPTION( "Method " #module "::" #function                             \
                                 " also not available on " << *loc );                         \
        }                                                                                     \
    }


#define LAMA_INTERFACE_FN_TT( function, loc, module, structname, T1, T2 )                     \
    typename module##Interface::structname<T1,T2>::function function;                         \
    function = loc->getInterface().module.function<T1,T2>();                                  \
    if ( function == NULL )                                                                   \
    {                                                                                         \
        LAMA_THROWEXCEPTION( "Method " #module "::" #function " not available on " << *loc ); \
    }

namespace lama
{

/** Interface for BLAS routines.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 *
 *  The interface is grouped in several structures 
 *
 *  - BLAS1 for BLAS level-1 routines
 *  - BLAS2 for BLAS level-2 routines
 *  - BLAS3 for BLAS level-3 routines
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

    // declare variables of function pointers, one for each type

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

    BLASInterface();
};

/** Structure with pointers for all Utils methods. */

struct UtilsInterface
{
    /** @brief Structure with type definitions for methods on indexes */

    struct Indexes
    {
        /** Check that all values in array are in a certain range.
         *
         *  @param array is an array of index values
         *  @param n is the size of array
         *  @param size specifies the range in which array values must fit
         *  @return true if \f$ 0 \le array[i] < size for all i = 0, ..., n-1\f$
         */

        typedef bool ( *validIndexes) ( const IndexType array[], const IndexType n, const IndexType size );
    };

    /** Structure with type defintions for all reduction methods.
     *
     *  @tparam ValueType specifies the value type used in the reduction.
     */

    template<typename ValueType>
    struct Reductions
    {
        /** Sum n contiguously stored values.
         *
         *  @param array is an array of values
         *  @param n is the size of array
         *  @return sum of all values in array
         */

        typedef ValueType ( *sum ) ( const ValueType array[], const IndexType n );

        typedef ValueType ( *maxval ) ( const ValueType array[], const IndexType n );

        /** Absolute maximum of n contiguously stored values */

        typedef ValueType ( *absMaxVal ) ( const ValueType array[], const IndexType n );

        /** Building absolute maximum of element-wise difference of vector elements.
         *
         *  @param array1i[in] first array
         *  @param array2i[in] second array
         *  @param n           size of array1 and array2
         *  @returns           max( abs( array1[i] - array2[i] ) ), \f$ 0 \le i < n \f$
         *
         *  Function is helpful to compute maximum norm for vectors and matrices
         */

        typedef ValueType ( *absMaxDiffVal ) ( const ValueType array1[], const ValueType array2[], const IndexType n );

        /** Predicate that tests whether a sequene is sorted.
         *
         *  @param[in] array values to be checked
         *  @param[in] n number of values to check
         *  @param[in] ascending if true check for ascending order, otherwise for descending
         */

        typedef bool ( *isSorted) ( const ValueType array[], const IndexType n, bool ascending );
    };

    template<typename ValueType>
    struct Setter
    {
        /** Set all elements of a contiguous array with a value. */

        typedef void ( *setVal ) ( ValueType array[], const IndexType n, const ValueType val );

        /** Set all elements of a contiguous array with its order number 0, 1, 2, ... */

        typedef void ( *setOrder ) ( ValueType array[], const IndexType n );
    };

    template<typename ValueType>
    struct Getter
    {
        typedef ValueType ( *getValue ) ( const ValueType* array, const IndexType i );
    };

    LAMA_INTERFACE_DEFINE_T( Getter, getValue )

    template<typename ValueType1, typename ValueType2>
    struct Copy
    {
        /** Set out[i] = in[i],  0 <= i < n */

        typedef void ( *set ) ( ValueType1 out[],
                                const ValueType2 in[],
                                const IndexType n );

        /** Set out[i] = in[ indexes[i] ],  \f$0 \le i < n\f$ */

        typedef void ( *setGather ) ( ValueType1 out[],
                                      const ValueType2 in[],
                                      const IndexType indexes[],
                                      const IndexType n );

        /** Set out[ indexes[i] ] = in [i] */

        typedef void ( *setScatter ) ( ValueType1 out[],
                                       const IndexType indexes[],
                                       const ValueType2 in[],
                                       const IndexType n );
    };

    template<typename ValueType>
    struct Math
    {
        /** Set array[i] = 1.0 / array[i],  0 <= i < n
         *
         *  @param[in,out] array is the array to invert
         *  @param         n     is the number of entries to invert
         */

        typedef void ( *invert ) ( ValueType array[], const IndexType n );
    };

    template<typename ValueType, typename OtherValueType>
    struct Transform
    {
        typedef void ( *scale ) ( ValueType values[],
                                  const IndexType n,
                                  const OtherValueType value );
    };

    LAMA_INTERFACE_DEFINE( Indexes, validIndexes )

    LAMA_INTERFACE_DEFINE_T( Reductions, sum )
    LAMA_INTERFACE_DEFINE_T( Reductions, maxval )
    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxVal )
    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxDiffVal )
    LAMA_INTERFACE_DEFINE_T( Reductions, isSorted )

    LAMA_INTERFACE_DEFINE_T( Setter, setVal )
    LAMA_INTERFACE_DEFINE_T( Setter, setOrder )

    LAMA_INTERFACE_DEFINE_TT( Copy, setGather )
    LAMA_INTERFACE_DEFINE_TT( Copy, setScatter )
    LAMA_INTERFACE_DEFINE_TT( Copy, set )

    LAMA_INTERFACE_DEFINE_T( Math, invert )

    LAMA_INTERFACE_DEFINE_TT( Transform, scale )

    /** Constructor initializes all function pointers with nullPtr */

    UtilsInterface ();
};

/** Interface for utility functions to be used in CSR storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct CSRUtilsInterface
{
    /** Structure with type definitions for spmv */

    template<typename ValueType>
    struct Operations
    {
        typedef void ( *sortRowElements ) ( IndexType csrJA[],
                                            ValueType csrValues[],
                                            const IndexType csrIA[],
                                            const IndexType numRows,
                                            const bool diagonalFlag ) ;
    };

    LAMA_INTERFACE_DEFINE_T( Operations, sortRowElements )

    /** Structure with type definitions for solver routines */

    template<typename ValueType>
    struct Solver
    {
       /** Method to compute one iteration step in Jacobi method
        *
        *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
        *
        */
       typedef void ( *jacobi ) ( ValueType* const solution,
                                  const IndexType csrIA[],
                                  const IndexType csrJA[],
                                  const ValueType csrValues[],
                                  const ValueType oldSolution[],
                                  const ValueType rhs[],
                                  const ValueType omega,
                                  const IndexType numRows,
                                  class SyncToken* syncToken );

        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         */
        typedef void ( *jacobiHalo ) ( ValueType solution[],
                                       const IndexType localIA[],
                                       const ValueType localValues[],
                                       const IndexType haloIA[],
                                       const IndexType haloJA[],
                                       const ValueType haloValues[],
                                       const IndexType haloRowIndexes[],
                                       const ValueType oldSolution[],
                                       const ValueType omega,
                                       const IndexType numNonEmptyRows );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )

    /** Structure with type definitions for offset routines. */

    struct Offsets
    {
        /** This method makes an offset array from the sizes.
         *
         *  @param[in,out] array contains counter values and later the offsets
         *  @param[in]    n is the number of values, array must contain one additional value
         *  @returns      the total number of values
         *
         *  \code
         *    array  :    3    7   8   4   2  x
         *    array  :    0   10  15  12  16  18  -> returns 18
         *  \endcode
         *
         *  Important: sizes must have numRows + 1 allocated entries.
         *
         */

        typedef IndexType ( *sizes2offsets ) ( IndexType array[], const IndexType n );

        /** This method computes size array from an offset array.
         *
         *  @param[out] sizes will contain the sizes (e.g. for each row ), has numRows entries
         *  @param[in] offsets contains the offsets, has numRows + 1 entries
         *  @param[in] n is size of array sizes, n + 1 is size of array offsets
         *
         *  \code
         *           offsets :   0   5  11   16  19  19  21  23
         *           sizes   :   5   6   5    3   0   2   2
         *  \endcode
         *
         */

        typedef void ( *offsets2sizes ) ( IndexType sizes[], const IndexType offsets[], const IndexType n );

        /** check for a legal offset array. */

        typedef bool ( *validOffsets ) ( const IndexType array[], const IndexType n, const IndexType total );

        /** This method computes the row sizes for result matrix C of matrix add A + B
         *
         *  @param[out] cSizes array of length numRows, will contain number of entries
         *  @param[in]  numRows number of rows for matrix A and B
         *  @param[in]  numColumns number of columns for matrix A and B
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         *
         *  Note: filling the result matrix must use the same flag for diagonalProperty
         *        otherwise the row sizes/offsets will not match
         */

        typedef IndexType ( *matrixAddSizes ) ( IndexType cIa[], const IndexType numRows,
                                           const IndexType numColumns, bool diagonalProperty,
                                           const IndexType aIA[], const IndexType aJA[],
                                           const IndexType bIA[], const IndexType bJA[] );

        /** This method computes the row sizes for result matrix C of matrix multiplication A x B
         *
         *  @param[out] cSizes array of length numRows, will contain number of entries
         *  @param[in]  numRows number of rows for matrix C and A
         *  @param[in]  numColumns number of columns for matrix C and B
         *  @param[in]  diagonalProperty if true, diagonal elements will count in any case
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         */

        typedef IndexType ( *matrixMultiplySizes ) ( IndexType cIa[], const IndexType numRows,
                                                     const IndexType numColumns, bool diagonalProperty,
                                                     const IndexType aIA[], const IndexType aJA[],
                                                     const IndexType bIA[], const IndexType bJA[] );

        /** This method computes the column indexes for result matrix C of matrix multiplication A x B
         *
         *  @param[out] cJA array will contain the column indexes, size is cIA[numRows]
         *  @param[in]  CIA array with row offsets as computed by matrixMultiplySizes + sizes2offsets
         *  @param[in]  numRows number of rows for matrix C and A
         *  @param[in]  numColumns number of columns for matrix C and B
         *  @param[in]  diagonalProperty if true, diagonal elements will filled in at begin of each row
         *  @param[in]  aIA, aJA are the index arrays of matrix A
         *  @param[in]  bIA, bJA are the index arrays of matrix B
         */

        typedef void ( *matrixMultiplyJA ) ( IndexType cJA[], const IndexType cIA[],
                                             const IndexType numRows, const IndexType numColumns,
                                             bool diagonalProperty,
                                             const IndexType aIA[], const IndexType aJA[],
                                             const IndexType bIA[], const IndexType bJA[] );

        typedef bool ( *hasDiagonalProperty ) ( const IndexType numDiagonals,
                                                const IndexType csrIA[],
                                                const IndexType csrJA[] );

    };

    LAMA_INTERFACE_DEFINE( Offsets, sizes2offsets )
    LAMA_INTERFACE_DEFINE( Offsets, offsets2sizes )
    LAMA_INTERFACE_DEFINE( Offsets, validOffsets )
    LAMA_INTERFACE_DEFINE( Offsets, matrixAddSizes )
    LAMA_INTERFACE_DEFINE( Offsets, matrixMultiplySizes )
    LAMA_INTERFACE_DEFINE( Offsets, matrixMultiplyJA )
    LAMA_INTERFACE_DEFINE( Offsets, hasDiagonalProperty )

    /** Define structure that contains type definitions for the function pointers.
     *
     *  @tparam ValueType specifies the value type used in the arrays.
     *
     *  The structure is needed as type definition templates are unsupported in C++.
     */

    template<typename ValueType>
    struct Transpose
    {
        /** Function pointer for CSR to CSC conversion routine.
         *
         *  @param[out] cscIA, cscJA, cscValues is CSC output data
         *  @param[in]  csrIA, csrJA, csrValues is CSR input data
         *  @param numRows x numColumns is shape of input CSR matrix
         *
         *  Arrays must be big enough for the corresponing sizes.
         */

        typedef void( *convertCSR2CSC ) ( IndexType cscIA[],
                                          IndexType cscJA[],
                                          ValueType cscValues[],
                                          const IndexType csrIA[],
                                          const IndexType csrJA[],
                                          const ValueType csrValues[],
                                          IndexType numRows, IndexType numColumns,
                                          IndexType numValues );
    };

    LAMA_INTERFACE_DEFINE_T( Transpose, convertCSR2CSC )

    /** Define structure for multiplication routines.  */

    template<typename ValueType1, typename ValueType2>
    struct Scale
    {

        /** This operation multiplies each row with an own value.
         *
         *  csr[i,j] *= values[i], for i = 0, ..., numRows-1
         */

        typedef void ( *scaleRows ) ( ValueType1 csrValues[],
                                      const IndexType csrIA[],
                                      const IndexType numRows,
                                      const ValueType2 values[] );
    };

    LAMA_INTERFACE_DEFINE_TT( Scale, scaleRows )

    /** Structure with type definitions for reduction routines. */

    template <typename ValueType>
    struct Reductions
    {
        /** Get the maximal element-wise difference for two CSR matrices.
         *
         *  @param[in] numRows is number of rows for both matrices
         *  @param[in] sortedRows if true column indexes in JA arrays are sorted
         *  @param[in] csrIA1, csrJA1, csrValues1 is storage data of first matrix
         *  @param[in] csrIA2, csrJA2, csrValues2 is storage data of second matrix
         *  @returns maximal value of absolute difference between two matrix elements
         */

        typedef ValueType ( *absMaxDiffVal ) ( IndexType numRows, bool sortedRows,
           const IndexType csrIA1[], const IndexType csrJA1[], const ValueType csrValues1[],
           const IndexType csrIA2[], const IndexType csrJA2[], const ValueType csrValues2[] );
    };

    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxDiffVal )

    /** Define structure for multiplication routines.  */

    template<typename ValueType>
    struct Mult
    {

        /** This operation multiplies each row with an own value.
         *
         *  csr[i,j] *= values[i], for i = 0, ..., numRows-1
         */

        typedef void ( *scaleRows ) ( ValueType csrValues[],
                                      const IndexType csrIA[],
                                      const IndexType numRows,
                                      const ValueType values[] );

        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const ValueType beta,
                                       const ValueType y[],
                                       const IndexType numRows,
                                       const IndexType csrIA[],
                                       const IndexType csrJA[],
                                       const ValueType csrValues[],
                                       class SyncToken* syncToken );

        /** result = alpha * CSR-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param csrIA, csrJA, csrValues are arrays of CSR storage
         *  @param syncToken optional, if available starts asynchronous computation
         *
         *  Note: this routine does not provide the term 'beta * y' as it would require
         *        to run over the full result vector
         */

        typedef void ( *sparseGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const IndexType numNonZeroRows,
                                       const IndexType rowIndexes[],
                                       const IndexType csrIA[],
                                       const IndexType csrJA[],
                                       const ValueType csrValues[],
                                       class SyncToken* syncToken );

        /**  This method computes result = alpha * CSR * x + beta * y  with dense result, x, y
         *
         *   param[out] result  has size m x n
         *   param[in]  alpha scaling factor
         *   param[in]  x has size p x n
         *   param[in]  beta scaling factor
         *   param[in]  y has size m x n
         *   param[in]  csrIA offset array of CSR matrix, has size m + 1
         *   param[in]  csrJA has size csrIA[m], values between 0 and p-1
         *   param[in]  csrVaues is value array of CSR matrix
         *   param[in,out] syncToken might be used for asynchronous execution
         */

        typedef void ( *gemm ) ( ValueType result[],
                                 const ValueType alpha,
                                 const ValueType x[],
                                 const ValueType beta,
                                 const ValueType y[],
                                 const IndexType m,
                                 const IndexType n,
                                 const IndexType p,
                                 const IndexType csrIA[],
                                 const IndexType csrJA[],
                                 const ValueType csrValues[],
                                 class SyncToken* syncToken );

        /** computes c = alpha * a + beta * b for CSR sparse matrices a, b, c
         *
         *  @param[out] cJA, cValues are the matrix values of output matrix
         *  @param[in]  cIA contains already computed offsets
         *  @param[in]  numRows is number of rows for matrices a, b, c
         *  @param[in]  numColums is number of columns for matrices a, b, c
         *  @param[in]  diagonalProperty if true result matrix gets diagonal property
         *  @param[in]  alpha is a scalar scaling factor for matrix a
         *  @param[in]  aIA, aJA, aValues is input matrix a in CSR format
         *  @param[in]  beta is a scalar scaling factor for matrix b
         *  @param[in]  bIA, bJA, bValues is input matrix b in CSR format
         *
         *  In this routine the row offsets of C must already be determined
         *  before.
         */

        typedef void ( *matrixAdd ) ( IndexType cJA[],
                                      ValueType cValues[],
                                      const IndexType cIA[],
                                      const IndexType numRows,
                                      const IndexType numColumns,
                                      const bool diagonalProperty,
                                      const ValueType alpha,
                                      const IndexType aIA[],
                                      const IndexType aJA[],
                                      const ValueType aValues[],
                                      const ValueType beta,
                                      const IndexType bIA[],
                                      const IndexType bJA[],
                                      const ValueType bValues[] );

        /** computes c = alpha * a * b for CSR sparse matrices a, b, c
         *
         *  @param[out] cValues are the matrix values of output matrix
         *  @param[in]  cIA, cJA contain structure of output matrix
         *  @param[in]  numRows is number of rows for matrix c and a
         *  @param[in]  alpha is a scalar scaling factor
         *  @param[in]  aIA, aJA, aValues is input matrix a in CSR format
         *  @param[in]  bIA, bJA, bValues is input matrix b in CSR format
         *
         *  In this routine the final structure of C must already be determined
         *  before. Only those values of C are computed for which cIA, cJA specify
         *  an available entry.
         */

        typedef void ( *matrixMultiply ) ( const IndexType cIa[],
                                           IndexType cJA[],
                                           ValueType cValues[],
                                           const IndexType numRows,
                                           const IndexType numColumns,
                                           const ValueType alpha,
                                           bool diagonalProperty,
                                           const IndexType aIA[],
                                           const IndexType aJA[],
                                           const ValueType aValues[],
                                           const IndexType bIA[],
                                           const IndexType bJA[],
                                           const ValueType bValues[] );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, gemm )
    LAMA_INTERFACE_DEFINE_T( Mult, sparseGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, matrixAdd )
    LAMA_INTERFACE_DEFINE_T( Mult, matrixMultiply )

    CSRUtilsInterface ();
};

/** Interface for utility functions to be used in Dense storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct DenseUtilsInterface
{
    /** Function pointer type definitions for dense storage. */

    template<typename DenseValueType>
    struct Counting
    {
         /** Counting non-zero values in dense storage for conversion to CSR
          *
          *  @param[out] csrSizes is an array that contains for each row the number of non-zero elements
          *  @param[in]  diagonalFlag if true the diagonal elements are counted in any case
          *  @param[in]  numRows number of rows
          *  @param[in]  numColumns number of columns
          *  @param[in]  denseValues size is numRows x numColumns, array with all matrix elements of dense format
          *  @param[in]  eps is threshold when an element is to be considered as non-zero
          *
          *  The matrix values are stored row-wise in denseValues.
          */

        typedef void ( *getCSRSizes )( IndexType csrSizes[],
                                       bool diagonalFlag,
                                       const IndexType numRows,
                                       const IndexType numColumns,
                                       const DenseValueType denseValues[],
                                       const DenseValueType eps );
    };

    LAMA_INTERFACE_DEFINE_T( Counting, getCSRSizes )

    /** Function pointer type definitions for conversion on dense storage. */

    template<typename DenseValueType, typename CSRValueType>
    struct Conversions
    {
        /** Convesion of dense matrix to CSR storage format
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in]  csrIA is the array with the offsets (must already be available before)
         *  @param[in]  diagonalFlag if true the diagonal elements are also filled
         *  @param[in]  numRows number of rows
         *  @param[in]  numColumns number of columns
         *  @param[in]  denseValues size is numRows x numColumns, array with all matrix elements of dense format
         *  @param[in]  eps is threshold when an element is to be considered as non-zero
         *
         *  Very important: the offsets in csrIA must correspond to the csrSizes computed
         *                  by getCSRSizes.
         */
        typedef void ( *getCSRValues ) ( IndexType csrJA[],
                                         CSRValueType csrValues[],
                                         const IndexType csrIA[],
                                         const bool diagonalFlag,
                                         const IndexType numRows,
                                         const IndexType numColumns,
                                         const DenseValueType denseValues[],
                                         const DenseValueType eps );

        /** Conversion of CSR format to dense matrix. */

        typedef void ( *setCSRValues ) ( DenseValueType denseValues[],
                                         const IndexType numRows,
                                         const IndexType numColumns,
                                         const IndexType csrIA[],
                                         const IndexType csrJA[],
                                         const CSRValueType csrValues[] );
    };

    LAMA_INTERFACE_DEFINE_TT( Conversions, setCSRValues )
    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )

    /** Function pointer type definitions for copying on dense storage. */

    template<typename DenseValueType1, typename DenseValueType2>
    struct Copy
    {
        /** Copy values of dense matrix; supports also conversion. */

        typedef void ( *copyDenseValues ) ( DenseValueType1 newValues[],
                                            const IndexType numRows,
                                           const IndexType numColumns,
                                           const DenseValueType2 oldValues[] );

         /** Get diagonal of a dense matrix, type conversion is supported. */

         typedef void ( *getDiagonal ) ( DenseValueType1 diagonalValues[],
                                         const IndexType numDiagonalValues,
                                         const DenseValueType2 denseValues[],
                                         const IndexType numRows,
                                         const IndexType numColumns );

         /** Set diagonal of a dense matrix, type conversion is supported. */

         typedef void ( *setDiagonal ) ( DenseValueType1 denseValues[],
                                         const IndexType numRows,
                                         const IndexType numColumns,
                                         const DenseValueType2 diagonalValues[],
                                         const IndexType numDiagonalValues );
    };

    LAMA_INTERFACE_DEFINE_TT( Copy, copyDenseValues )
    LAMA_INTERFACE_DEFINE_TT( Copy, setDiagonal )
    LAMA_INTERFACE_DEFINE_TT( Copy, getDiagonal )

    /** Function pointer type definitions for modification of dense storage. */

    template<typename DenseValueType>
    struct Modify
    {
        /** Scale all elements of the dense matrix with a value */

        typedef void ( *scaleValue ) ( DenseValueType denseValues[],
                                       const IndexType numRows,
                                       const IndexType numColumns,
                                       const DenseValueType val );

        /** Set diagonal elements with one and the same value. */

        typedef void ( *setDiagonalValue ) ( DenseValueType denseValues[],
                                             const IndexType numRows,
                                             const IndexType numColumns,
                                             const DenseValueType val );
    };

    LAMA_INTERFACE_DEFINE_T( Modify, scaleValue )
    LAMA_INTERFACE_DEFINE_T( Modify, setDiagonalValue )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numColumns is the number of columns, size of x
         *  @param numValues, array with the dense matrix values
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const ValueType beta,
                                       const ValueType y[],
                                       const IndexType numRows,
                                       const IndexType numColumns,
                                       const ValueType denseValues[],
                                       class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )

    /** Default constructor initializes all function pointers with NULL. */

    DenseUtilsInterface ();
};

/** Interface for utility functions to be used in ELL storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct ELLUtilsInterface
{
    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct Solver
    {
        /**
         *
         *  @param[out] solution is solution vector, size is numRows
         *  @param[in]  numRows is size of vectors and number of rows for matrix
         *  @param[in]  ellNumValuesPerRow is maximal number of non-zero entries
         *  @param[in]  ellSizes, ellJA, ellValues are arrays of ELL storage, numRows x numRows
         *  @param[in]  oldSolution is the old solution, size is numRows
         *  @param[in]  rhs is right hand side vector, size is numRows
         *  @param[in]  omega is scaling factor
         *  @param[in]  syncToken optional, NULL synchronous execution, otherwise asynchronous
         *
         *  The ELL storage stands for a square matrix and must have diagonal property.
         */

        typedef void ( *jacobi ) ( ValueType solution[],
                                   const IndexType numRows,
                                   const IndexType ellNumValuesPerRow,
                                   const IndexType ellSizes[],
                                   const IndexType ellJA[],
                                   const ValueType ellValues[],
                                   const ValueType oldSolution[],
                                   const ValueType rhs[],
                                   const ValueType omega,
                                   class SyncToken* syncToken );

        typedef void ( *jacobiHalo ) (ValueType solution[],
                                      const IndexType numRows,
                                      const ValueType diagonal[],
                                      const IndexType ellNumValuesPerRow,
                                      const IndexType ellSizes[],
                                      const IndexType ellJA[],
                                      const ValueType ellValues[],
                                      const IndexType rowIndexes[],
                                      const IndexType numNonEmptyRows,
                                      const ValueType oldSolution[],
                                      const ValueType omega,
                                      class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )

    /** Conversion routines between ELL and CSR storage format. */

    template<typename ELLValueType, typename CSRValueType>
    struct Conversions
    {
        /** Conversion from ELL data to CSR data
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in]  csrIA is the array with the offsets (must already be available before)
         *  @param[in]  numRows is number of rows
         *  @param[in]  ellSizes is the number of values in each row
         *  @param[in]  ellJA are the column indexes for ELL format
         *  @param[in]  ellValues are the stored matrix values for ELL format
         */

        typedef void ( *getCSRValues ) ( IndexType csrJA[],
                                         CSRValueType csrValues[],
                                         const IndexType csrIA[],
                                         const IndexType numRows,
                                         const IndexType ellSizes[],
                                         const IndexType ellJA[],
                                         const ELLValueType ellValues[] );

        /** Conversion from CSR data to ELL data      */

        typedef void( *setCSRValues ) ( IndexType ellJA[],
                                        ELLValueType ellValues[],
                                        const IndexType ellSizes[],
                                        const IndexType numRows,
                                        const IndexType numValuesPerRow,
                                        const IndexType csrIA[],
                                        const IndexType csrJA[],
                                        const CSRValueType csrValues[] );

    };

    LAMA_INTERFACE_DEFINE_TT( Conversions, setCSRValues )
    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )

    template<typename ValueType>
    struct Helper
    {

        /** Compresses the given IA array using the values array and epsilon
         *
         * @param[in]  IA that should be compressed
         * @param[in]  related JA array
         * @param[in]  related values array
         * @param[in]  number of rows
         * @param[in]  epsilon
         * @param[out] new created IA
         */

        typedef void ( *compressIA ) ( const IndexType IA[],
                                       const IndexType JA[],
                                       const ValueType values[],
                                       const IndexType numRows,
                                       const ValueType eps,
                                       IndexType newIA[] );

        /** Compresses the given JA and values array using epsilon
         *
         * @param[in]  IA that should be compressed
         * @param[in]  related JA array
         * @param[in]  related values array
         * @param[in]  number of rows
         * @param[in]  epsilon
         * @param[out] new created JA
         * @param[out] new created values
         */

        typedef void ( *compressValues ) ( const IndexType IA[],
                                           const IndexType JA[],
                                           const ValueType values[],
                                           const IndexType numRows,
                                           const ValueType eps,
                                           IndexType newJA[],
                                           ValueType newValues[] );

    };

    LAMA_INTERFACE_DEFINE_T( Helper, compressIA )
    LAMA_INTERFACE_DEFINE_T( Helper, compressValues )

    template<typename ValueType, typename OtherValueType>
    struct Getter
    {
        /** Returns a row of the matrix as dense vector
         *
         *  @param[out] row as dense vector that will be returned
         *  @param[in]  i is the row that should be returned
         *  @param[in]  numRows is the number of rows of the ELL matrix
         *  @param[in]  numColums is size of ia
         *  @param[in]  ia is the ELL sizes array
         *  @param[in]  ja is the ELL ja array
         *  @param[in]  values is the ELL values array
         */

        typedef void ( *getRow ) ( OtherValueType *row,
                                   const IndexType i,
                                   const IndexType numRows,
                                   const IndexType numColumns,
                                   const IndexType *ia,
                                   const IndexType *ja,
                                   const ValueType *values );

        /** Returns one element of the matrix
         *
         *  @param[in] i is the row of the returned element
         *  @param[in] j is the column of the returned element
         *  @param[in] numRows is the number of rows of the matrix
         *  @param[in] ia is the ELL sizes array
         *  @param[in] ja is the ELL ja array
         *  @param[in] values is the ELL values array
         */

        typedef OtherValueType ( *getValue ) ( const IndexType i,
                                          const IndexType j,
                                          const IndexType numRows,
                                          const IndexType *ia,
                                          const IndexType *ja,
                                          const ValueType *values );
    };

    LAMA_INTERFACE_DEFINE_TT( Getter , getRow )
    LAMA_INTERFACE_DEFINE_TT( Getter, getValue )

    struct Operations
    {
        typedef IndexType ( *countNonEmptyRowsBySizes ) ( const IndexType sizes[],
                                                          const IndexType numRows );

        typedef void ( *setNonEmptyRowsBySizes ) (  IndexType rowIndexes[],
                                                    const IndexType numNonEmptyRows,
                                                    const IndexType sizes[],
                                                    const IndexType numRows );

        typedef bool ( *hasDiagonalProperty ) ( const IndexType numDiagonals,
                                                        const IndexType ellJA[] );

        typedef void ( *check ) ( const IndexType mNumRows,
                                  const IndexType mNumValuesPerRow,
                                  const IndexType mNumColumns,
                                  const IndexType *ia,
                                  const IndexType *ja,
                                  const char* msg );
    };

    LAMA_INTERFACE_DEFINE( Operations, countNonEmptyRowsBySizes )
    LAMA_INTERFACE_DEFINE( Operations, setNonEmptyRowsBySizes )
    LAMA_INTERFACE_DEFINE( Operations, hasDiagonalProperty )
    LAMA_INTERFACE_DEFINE( Operations, check )

    /** Define structure for multiplication routines.  */

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * ELL-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param ellIA, ellJA, csrValues are arrays of ELL storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const ValueType beta,
                                       const ValueType y[],
                                       const IndexType numRows,
                                       const IndexType numNonZerosPerRows,
                                       const IndexType ellIA[],
                                       const IndexType ellJA[],
                                       const ValueType ellValues[],
                                       class SyncToken* syncToken );

        /** result = alpha * ELL-Matrix * x, CSR matrix has only some non-zero rows
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param numNonZeroRows is size of rowIndexes
         *  @param rowIndexes are indexes of non-empty rows in matrix
         *  @param ellIA, ellJA, csrValues are arrays of ELL storage
         *  @param syncToken optional, if available starts asynchronous computation
         *
         *  Note: this routine does not provide the term 'beta * y' as it would require
         *        to run over the full result vector
         */

        typedef void ( *sparseGEMV ) ( ValueType result[],
                                       const IndexType numRows,
                                       const IndexType numNonZerosPerRows,
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const IndexType numNonZeroRows,
                                       const IndexType rowIndexes[],
                                       const IndexType ellIA[],
                                       const IndexType ellJA[],
                                       const ValueType ellValues[],
                                       class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )
    LAMA_INTERFACE_DEFINE_T( Mult, sparseGEMV )

    /** Structure with type definitions for reduction routines */

    template<typename ValueType>
    struct Reductions
    {
         /** This method returns the maximal absolute value of an ELLPACK matrix. */

         typedef ValueType ( *absMaxVal ) ( const IndexType numRows,
                                            const IndexType numNonZerosPerRow,
                                            const IndexType ellIA[],
                                            const ValueType ellValues[]
                                          );
    };

    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxVal )

    template<typename ValueType, typename OtherValueType>
    struct Scale
    {
        typedef void ( *scaleValue ) ( const IndexType numRows,
                                       const IndexType Ia[],
                                       ValueType mValues[],
                                       const OtherValueType values[] );
    };

    LAMA_INTERFACE_DEFINE_TT( Scale, scaleValue )

    template<typename ValueType>
    struct MatrixTimesMatrix
    {
        /** Method to compute the resulting IA array of an matrix times matrix multiplication
         *
         *  @param IA array of the left hand matrix
         *  @param JA array of the left hand matrix
         *  @param number of rows of the left hand matrix
         *  @param IA array of the right hand matrix
         *  @param JA array of the right hand matrix
         *  @param number of rows of the right hand matrix
         *  @param resulting IA array
         *
         */

        typedef void ( *computeIA ) ( const IndexType aIA[],
                                      const IndexType aJA[],
                                      const IndexType aNumRows,
                                      const IndexType bIA[],
                                      const IndexType bJA[],
                                      const IndexType bNumRows,
                                      IndexType cIA[] );

        /** Method to compute the resulting JA and Values arrays of the matrix times matrix multiplication
         *
         * @param IA array of the left hand matrix
         * @param JA array of the left hand matrix
         * @param values array of the left hand matrix
         * @param number of rows of the left hand matrix
         * @param IA array of the right hand matrix
         * @param JA array of the right hand matrix
         * @param values array of the right hand matrix
         * @param number of rows of the right hand matrix
         * @param alpha of the multiplication
         * @param resulting JA array
         * @param resulting values array
         */

        typedef void ( *computeValues ) ( const IndexType aIA[],
                                          const IndexType aJA[],
                                          const ValueType aValues[],
                                          const IndexType aNumRows,
                                          const IndexType bIA[],
                                          const IndexType bJA[],
                                          const ValueType bValues[],
                                          const IndexType bNumRows,
                                          const ValueType alpha,
                                          const IndexType cIA[],
                                          IndexType cJA[],
                                          ValueType cValues[] );

        /** Method to compute the resulting IA array of an matrix adding
         *
         *  @param IA array of the left hand matrix
         *  @param JA array of the left hand matrix
         *  @param number of rows of the left hand matrix
         *  @param IA array of the right hand matrix
         *  @param JA array of the right hand matrix
         *  @param number of rows of the right hand matrix
         *  @param resulting IA array
         *
         */

        typedef void ( *addComputeIA ) ( const IndexType aIA[],
                                         const IndexType aJA[],
                                         const IndexType aNumRows,
                                         const IndexType bIA[],
                                         const IndexType bJA[],
                                         const IndexType bNumRows,
                                         IndexType cIA[] );

        /** Method to compute the resulting JA and Values arrays of the matrix adding
         *
         * @param IA array of the left hand matrix
         * @param JA array of the left hand matrix
         * @param values array of the left hand matrix
         * @param number of rows of the left hand matrix
         * @param IA array of the right hand matrix
         * @param JA array of the right hand matrix
         * @param values array of the right hand matrix
         * @param number of rows of the right hand matrix
         * @param beta of the adding
         * @param resulting JA array
         * @param resulting values array
         */

        typedef void ( *addComputeValues ) ( const IndexType aIA[],
                                             const IndexType aJA[],
                                             const ValueType aValues[],
                                             const IndexType aNumRows,
                                             const IndexType bIA[],
                                             const IndexType bJA[],
                                             const ValueType bValues[],
                                             const IndexType bNumRows,
                                             const ValueType beta,
                                             const IndexType cIA[],
                                             IndexType cJA[],
                                             ValueType cValues[] );

    };

    LAMA_INTERFACE_DEFINE_T( MatrixTimesMatrix, computeIA )
    LAMA_INTERFACE_DEFINE_T( MatrixTimesMatrix, computeValues )
    LAMA_INTERFACE_DEFINE_T( MatrixTimesMatrix, addComputeIA )
    LAMA_INTERFACE_DEFINE_T( MatrixTimesMatrix, addComputeValues )

    ELLUtilsInterface ();

};

/** Interface for utility functions to be used in JDS storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct JDSUtilsInterface
{
    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct Solver
    {
       /** Method to compute one iteration step in Jacobi method
        *
        *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
        *
        */
       typedef void ( *jacobi ) ( ValueType* const solution,
                                  const IndexType numRows,
                                  const IndexType jdsPerm[],
                                  const IndexType jdsIlg[],
                                  const IndexType jdsNumDiagonals,
                                  const IndexType jdsDlg[],
                                  const IndexType jdsJA[],
                                  const ValueType jdsValues[],
                                  const ValueType oldSolution[],
                                  const ValueType rhs[],
                                  const ValueType omega,
                                  class SyncToken* syncToken );

        /** Method to compute one iteration step in Jacobi method with halo.  */

        typedef void ( *jacobiHalo ) ( ValueType solution[],
                                       const IndexType numRows,
                                       const ValueType invDiagonal[],
                                       const IndexType numDiagonals,
                                       const IndexType jdsHaloPerm[],
                                       const IndexType jdsHaloIlg[],
                                       const IndexType jdsHaloDlg[],
                                       const IndexType jdsHaloJA[],
                                       const ValueType jdsHaloValues[],
                                       const ValueType oldSolution[],
                                       const ValueType omega,
                                       class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )

    struct Sort
    {
        /** Stable sorting of values in array in descending order.
         *
         *  @param[in,out] array are the values to be sorted
         *  @param[in,out] perm, where perm[i] has the value of the original position
         *  @param[in]    n is the number of values to be sorted
         *
         *  \code
         *           array =   1  4   1  8  5  7
         *           perm  =   0  1   2  3  4  5
         +
         *           array =   8  7   5  4  1  1
         *           perm  =   3  5   4  1  0  2
         *  \endcode
         */

        typedef void ( *sortRows ) ( IndexType array[],
                                     IndexType perm[],
                                     const IndexType n );

        /** Compute the inverse permutation for a given permutation.
         *
         *  inversePerm [ perm [i] ] == i , 0 <= i < n
         *
         *  @param[out] inversePerm, size = n, will contain the inverse permutation
         *  @param[in] perm, size = n, is input permuation of 0, ..., n-1
         *  @param[in] n specifies the size of perm and inversePerm
         *
         *  /code
         *       perm      2  5  1  4  6  3  0
         *     inperm      6  2  0  5  3  1  4
         *  /endcode
         */

        typedef void ( *setInversePerm ) ( IndexType inversePerm[],
                                           const IndexType perm[],
                                           const IndexType n );

        /** Compute dlg array from ilg array.
         *
         *  @param[out] dlg is the array with sizes of the columns
         *  @param[in]  numDiagonals is size of dlg
         *  @param[in]  ilg is the array with sizes of the rows
         *  @param[in]  numRows is the number of rows, size of ilg
         *
         *  The values in ilg must be descending. The same will be true
         *  for the output array dlg.
         *
         *  /code
         *       ilg       4  3  2  2  1   dlg
         *       5         x  x  x  x  x
         *       4         x  x  x  x
         *       2         x  x
         *       1         x
         *  /endcode
         */

         typedef IndexType ( *ilg2dlg ) ( IndexType dlg[], const IndexType numDiagonals,
                                          const IndexType ilg[], const IndexType numRows );

    };

    LAMA_INTERFACE_DEFINE( Sort, sortRows )
    LAMA_INTERFACE_DEFINE( Sort, setInversePerm )
    LAMA_INTERFACE_DEFINE( Sort, ilg2dlg )

    template<typename JDSValueType, typename CSRValueType>
    struct Conversions
    {
        /** Conversion of JDS storage data to CSR data
         *
         *  @param[out]  csrJA is array with column indexes
         *  @param[out]  csrValues is array with non-zero values
         *  @param[in]   csrIA is offset array (must be computed before)
         *  @param[in]   numRows number of rows in matrix
         *  @param[in]   jdsPerm with jdsPerm[ii] is original index of row i
         *  @param[in]   jdsILG with size of entries in row i
         *  @param[in]   jdsDLG distances of columns
         *  @param[in]   jdsJA column indexes
         *  @param[in]   jdsValues matrix values
         */

        typedef void ( *getCSRValues ) ( IndexType csrJA[],
                                         CSRValueType csrValues[],
                                         const IndexType csrIA[],
                                         const IndexType numRows,
                                         const IndexType jdsPerm[],
                                         const IndexType jdsILG[],
                                         const IndexType jdsDLG[],
                                         const IndexType jdsJA[],
                                         const JDSValueType jdsValues[] );

        /** Conversion of CSR storage data to JDS data
         *
         *  @param[out]  jdsJA column indexes
         *  @param[out]  jdsValues matrix values
         *  @param[in]   numRows number of rows in matrix
         *  @param[in]   jdsPerm with jdsPerm[ii] is original index of row i
         *  @param[in]   jdsILG with size of entries in row i
         *  @param[in]   jdsDLG distances of columns
         *  @param[in]   csrIA is offset array (must be computed before)
         *  @param[in]   csrJA is array with column indexes
         *  @param[in]   csrValues is array with non-zero values
         */

        typedef void( *setCSRValues ) ( IndexType jdsJA[],
                                        JDSValueType jdsValues[],
                                        const IndexType numRows,
                                        const IndexType jdsPerm[],
                                        const IndexType jdsILG[],
                                        const IndexType jdsDLG[],
                                        const IndexType csrIA[],
                                        const IndexType csrJA[],
                                        const CSRValueType csrValues[] );
    };

    // Define tables ( indexed by template value types) for all methods

    LAMA_INTERFACE_DEFINE_TT( Conversions, setCSRValues )
    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numDiagonals are number of diagonals, is size of jdsDLG
         *  @param jdsPerm, jdsILG, jdsDLG, jdsJA, jdsValues are arrays of JDS storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const ValueType beta,
                                       const ValueType y[],
                                       const IndexType numRows,
                                       const IndexType jdsPerm[],
                                       const IndexType jdsILG[],
                                       const IndexType ndlg,
                                       const IndexType jdsDLG[],
                                       const IndexType jdsJA[],
                                       const ValueType jdsValues[],
                                       class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )

    template<typename ValueType>
    struct Operations
    {
        typedef void ( *setDiagonalWithScalar ) ( const IndexType numDiagonal,
                                                  ValueType values[],
                                                  Scalar scalar );
    };

    LAMA_INTERFACE_DEFINE_T( Operations, setDiagonalWithScalar )

    template<typename ValueType, typename OtherValueType>
    struct Getter
    {
        typedef void ( *getRow ) ( OtherValueType row[],
                                   const IndexType i,
                                   const IndexType numColumns,
                                   const IndexType numRows,
                                   const IndexType perm[],
                                   const IndexType ilg[],
                                   const IndexType dlg[],
                                   const IndexType ja[],
                                   const ValueType values[] );

        typedef ValueType (*getValue ) ( const IndexType i,
                                         const IndexType j,
                                         const IndexType numRows,
                                         const IndexType* dlg,
                                         const IndexType* ilg,
                                         const IndexType* perm,
                                         const IndexType* ja,
                                         const ValueType* values );

    };

    LAMA_INTERFACE_DEFINE_TT( Getter, getRow )
    LAMA_INTERFACE_DEFINE_TT( Getter, getValue )

    template<typename ValueType, typename OtherValueType>
    struct Scale
    {
        typedef void ( *scaleValue ) ( const IndexType numRows,
                                       const IndexType perm[],
                                       const IndexType ilg[],
                                       const IndexType dlg[],
                                       ValueType mValues[],
                                       const OtherValueType values[] );
    };

    LAMA_INTERFACE_DEFINE_TT( Scale, scaleValue )

    struct Helper
    {
        typedef bool ( *checkDiagonalProperty ) ( const IndexType numDiagonals,
                                                  const IndexType numRows,
                                                  const IndexType numColumns,
                                                  const IndexType perm[],
                                                  const IndexType ja[],
                                                  const IndexType dlg[] );

        typedef bool ( *check ) ( const IndexType numRows,
                                  const IndexType numValues,
                                  const IndexType numColumns,
                                  const IndexType ja[],
                                  const IndexType ilg[],
                                  const IndexType dlg[] );
    };

    LAMA_INTERFACE_DEFINE( Helper, checkDiagonalProperty )
    LAMA_INTERFACE_DEFINE( Helper, check )

    JDSUtilsInterface ();
};

/** Interface for utility functions to be used in DIA storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct DIAUtilsInterface
{
    template<typename ValueType>
    struct Counting
    {
        /** Type definition of function pointer for counting sparse values in DIA storage
         *
         *  @param[out] csrSizes array with number of non-zero entries in each row
         *  @param[in]  diagonalFlag if true count also zero diagonal elements
         *  @param[in] numRows is the number of rows
         *  @param[in] numColumns is the number of columns
         *  @param[in] numDiagonals number of diagonals used in the DIA format
         *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
         *  @param[in] diaValues are stored values of the diagonals
         *  @param[in] eps threshold value when an element should be considered as zero
         *
         *  Note: the diagonals might contain zero entries so the number of non-zero
         *        elements might be less than number of stored elements
         *
         *  - csrSizes must have been allocated with at least numRows entries
         *  - diaOffsets has at least numDiagonals entries
         *  - diaValues has numDiagonals x max(numRows, numColumns) entries
         */

        typedef void ( *getCSRSizes ) ( IndexType csrSizes[],
                                        bool diagonalFlag,
                                        const IndexType numRows,
                                        const IndexType numColumns,
                                        const IndexType numDiagonals,
                                        const IndexType diaOffsets[],
                                        const ValueType diaValues[],
                                        const ValueType eps );

    };

    LAMA_INTERFACE_DEFINE_T( Counting, getCSRSizes )

    template<typename DIAValueType, typename CSRValueType>
    struct Conversions
    {
        /** Type definition of function pointer for conversion of DIA storage data to CSR data.
        *
        *  @param[out] csrJA will contain the column indexes
        *  @param[out] csrValues will contain the matrix elements
        *  @param[in] csrIA is the array with the offsets (must already be available before)
        *  @param[in] numRows is the number of rows
        *  @param[in] numColumns is the number of columns
        *  @param[in] numDiagonals number of diagonals used in the DIA format
        *  @param[in] diaOffsets diagonal offsets, size is numDiagonals
        *  @param[in] diaValues are stored values of the diagonals
        *  @param[in] eps threshold value when an element should be considered as zero
        *
        *   - csrIA has numRows + 1 entries
        *   - csrJA and csrValues must have at least numValues entries, numValues = csrIA[numRows]
        */

        typedef void ( *getCSRValues ) ( IndexType csrJA[],
                                         CSRValueType csrValues[],
                                         const IndexType csrIA[],
                                         const bool diagonalFlag,
                                         const IndexType numRows,
                                         const IndexType numColumns,
                                         const IndexType numDiagonals,
                                         const IndexType diaOffsets[],
                                         const DIAValueType diaValues[],
                                         const DIAValueType eps );
    };

    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param numValues is the number of diagonals in DIA storage
         *  @param diaOffsets, diaValues are arrays of DIA storage
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const ValueType beta,
                                       const ValueType y[],
                                       const IndexType numRows,
                                       const IndexType numColumns,
                                       const IndexType numDiagonals,
                                       const IndexType diaOffsets[],
                                       const ValueType diaValues[],
                                       class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )

    /** Structure with type definitions for solver routines */

    template<typename ValueType>
    struct Solver
    {
       /** Method to compute one iteration step in Jacobi method
        *
        *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
        *
        */
       typedef void ( *jacobi ) ( ValueType* const solution,
                                  const IndexType numColumns,
                                  const IndexType numDiagonals,
                                  const IndexType diaOffset[],
                                  const ValueType diaValues[],
                                  const ValueType oldSolution[],
                                  const ValueType rhs[],
                                  const ValueType omega,
                                  const IndexType numRows,
                                  class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )

    /** Structure with type definitions for reduction routines */

    template<typename ValueType>
    struct Reductions
    {
         /** This method returns the maximal absolute value of an ELLPACK matrix. */

         typedef ValueType ( *absMaxVal ) ( const IndexType numRows,
                                            const IndexType numColumns,
                                            const IndexType numDiagonals,
                                            const IndexType diaOffsets[],
                                            const ValueType diaValues[]
                                          );
    };

    LAMA_INTERFACE_DEFINE_T( Reductions, absMaxVal )

    /** Constructur of interface sets all function pointers to nullPtr. */

    DIAUtilsInterface ();
};

/** Interface for utility functions to be used in COO storage.
 *
 *  This interface contains function pointer type definitions for all used routines
 *  and tables with actual values for the functions.
 */

struct COOUtilsInterface
{
    struct Counting
    {
        /** Helper routine for conversion of COO format to CSR format to get sparse row sizes.
         *
         *  @param[out] csrSizes array with number of non-zero entries in each row
         *  @param[in] numRows number of rows
         *  @param[in] numValues number of non-zero values
         *  @param[in] array with row indexes of COO storage (size is numValues)
         */

        typedef void ( *getCSRSizes )( IndexType csrSizes[],
                                       const IndexType numRows,
                                       const IndexType numValues,
                                       const IndexType cooIA[] );
    };

    LAMA_INTERFACE_DEFINE( Counting, getCSRSizes )

    template<typename COOValueType, typename CSRValueType>
    struct Conversions
    {
        /** Helper routine for conversion COO to CSR
         *
         *  @param[out] csrJA will contain the column indexes
         *  @param[out] csrValues will contain the matrix elements
         *  @param[in] csrIA is the array with the offsets (must already be available before)
         *
         *   - csrIA has numRows + 1 entries
         *   - csrJA and csrValues must have at least numValues entries, numValues = csrIA[numRows]
         *
         *  Note: this routine preserves the diagonal property of the COO format
         */

        typedef void ( *getCSRValues )( IndexType csrJA[],
                                        CSRValueType csrValues[],
                                        IndexType csrIA[],
                                        const IndexType numRow,
                                        const IndexType numValues,
                                        const IndexType cooIA[],
                                        const IndexType cooJA[],
                                        const COOValueType cooValues[] );

        /** Conversion of CSR data to COO data.
         *
         *  @param[in] numRows is number of rows
         *  @param[in] numDiagonals is number of diagonal elements to be stored at the beginning
         *  @param[in] csrIA is CSR offset array
         *  @param[in] csrJA is CSR column indexes, size is csrIA[numRows]
         *  @param[in] csrValues is CSR non-zero elements, size is csrIA[numRows]
         *  @param[in] csrDiagonalProperty is true if CSR data has diagonal property
         *
         *  Be careful if COO data should have diagonal property, but CSR data has not
         */

        typedef void ( *setCSRValues ) ( IndexType cooIA[],
                                         IndexType cooJA[],
                                         COOValueType cooValues[],
                                         const IndexType numRows,
                                         const IndexType numDiagonals,
                                         const IndexType csrIA[],
                                         const IndexType csrJA[],
                                         const CSRValueType csrValues[],
                                         const bool csrDiagonalProperty );
    };

    LAMA_INTERFACE_DEFINE_TT( Conversions, getCSRValues )
    LAMA_INTERFACE_DEFINE_TT( Conversions, setCSRValues )

    template<typename ValueType>
    struct Mult
    {
        /** result = alpha * CSR-Matrix * x + b * y.
         *
         *  @param result is the result vector
         *  @param alpha is scaling factor for matrix x vector
         *  @param x is input vector for matrix multiplication
         *  @param beta is scaling factor for additional vector
         *  @param y is additional input vector to add
         *  @param numRows is number of elements for all vectors and rows of matrix
         *  @param cooIA, cooJA, cooValues are arrays of COO storage
         *  @param numValues is the size of the coo arrays
         *  @param syncToken optional, if available starts asynchronous computation
         */

        typedef void ( *normalGEMV ) ( ValueType result[],
                                       const ValueType alpha,
                                       const ValueType x[],
                                       const ValueType beta,
                                       const ValueType y[],
                                       const IndexType numRows,
                                       const IndexType cooIA[],
                                       const IndexType cooJA[],
                                       const ValueType cooValues[],
                                       const IndexType numValues,
                                       class SyncToken* syncToken );
    };

    LAMA_INTERFACE_DEFINE_T( Mult, normalGEMV )

    /** Structure with type definitions for solver routines */
    template<typename ValueType>
    struct Solver
    {
       /** Method to compute one iteration step in Jacobi method
        *
        *  solution = omega * ( rhs + B * oldSolution) * dinv  + ( 1 - omega ) * oldSolution
        *
        */
       typedef void ( *jacobi ) ( ValueType* const solution,
                                  const IndexType cooNumValues,
                                  const IndexType cooIA[],
                                  const IndexType cooJA[],
                                  const ValueType cooValues[],
                                  const ValueType oldSolution[],
                                  const ValueType rhs[],
                                  const ValueType omega,
                                  const IndexType numRows,
                                  class SyncToken* syncToken );

        /** Method to compute one iteration step in Jacobi method
         *
         *  solution -= omega * ( B(halo) * oldSolution) * dinv
         *
         */
        typedef void ( *jacobiHalo ) ( ValueType solution[],
                                       const ValueType diaValues[],
                                       const IndexType haloIA[],
                                       const IndexType haloJA[],
                                       const ValueType haloValues[],
                                       const IndexType haloRowIndexes[],
                                       const ValueType oldSolution[],
                                       const ValueType omega,
                                       const IndexType numNonEmptyRows );
    };

    LAMA_INTERFACE_DEFINE_T( Solver, jacobi )
    LAMA_INTERFACE_DEFINE_T( Solver, jacobiHalo )

    /** Default constructor initializes all function pointers with NULL. */

    COOUtilsInterface ();
};

/**
 * @brief LAMAInterface is the abstract class for Matrix-Vector operations.
 *
 * This class specifies all routines that have to be implemented on
 * a new device.
 *
 * All routines of this class are virtual and should be overwritten by derived
 * classes. The implementation of the methods of this base class will all throw
 * an exception.
 */
class LAMA_DLL_IMPORTEXPORT LAMAInterface : public Printable
{
public:

    LAMAInterface();

    virtual ~LAMAInterface();

    /** This method writes the name of the interface into the output stream.
     *  Output message helps to identify which interfaces have been used.
     *
     * @see Printable for more details.
     */
    virtual void writeAt( std::ostream& stream ) const;

    // one member variable for each Interface

    BLASInterface       BLAS;        //!< interface table for BLAS routines
    CSRUtilsInterface   CSRUtils;
    DenseUtilsInterface DenseUtils;
    ELLUtilsInterface   ELLUtils;
    JDSUtilsInterface   JDSUtils;
    DIAUtilsInterface   DIAUtils;
    COOUtilsInterface   COOUtils;
    UtilsInterface      Utils;

protected:

    LAMA_LOG_DECL_STATIC_LOGGER(logger)
};

} //namespace lama

#endif // LAMA_LAMA_INTERFACE_HPP_
