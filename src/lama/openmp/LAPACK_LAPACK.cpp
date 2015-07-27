/**
 * @file LAPACK_LAPACK.cpp
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
 * @brief LAPACK_LAPACK.cpp
 * @author lschubert
 * @date 02.07.2012
 * @since 1.0.0
 */

// hpp
#include <lama/openmp/LAPACK_LAPACK.hpp>

// others
#include <lama/openmp/BLAS_BLAS1.hpp>

#include <lama/LAMAInterfaceRegistry.hpp>
#include <lama/BLASInterface.hpp>
#include <tracing/tracing.hpp>
#include <common/unique_ptr.hpp>

// macros
#include <lama/macros/unused.hpp>

#include <lama/cblas.hpp>

//fallback if nothing is set in cmake
#if !defined(LAMA_FORTRAN_BLAS_STYLE_UNDERSCORE)
#if !defined(LAMA_FORTRAN_BLAS_STYLE_UPCASE)
#if !defined(LAMA_FORTRAN_BLAS_STYLE_LOWCASE)
#define LAMA_FORTRAN_BLAS_STYLE_LOWCASE
#endif
#endif
#endif

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#if defined( LAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )

#define F77_sgetrf sgetrf_
#define F77_dgetrf dgetrf_
#define F77_sgetri sgetri_
#define F77_dgetri dgetri_
#define F77_strtrs strtrs_
#define F77_dtrtrs dtrtrs_
#define F77_stptrs stptrs_
#define F77_dtptrs dtptrs_
#define F77_slaswp slaswp_
#define F77_dlaswp dlaswp_

#elif defined(LAMA_FORTRAN_BLAS_STYLE_UPCASE)

#define F77_sgetrf SGETRF
#define F77_dgetrf DGETRF
#define F77_sgetri SGETRI
#define F77_dgetri DGETRI
#define F77_strtrs STRTRS
#define F77_dtrtrs DTRTRS
#define F77_stptrs STPTRS
#define F77_dtptrs DTPTRS
#define F77_slaswp SLASWP
#define F77_dlaswp DLASWP

#elif defined(LAMA_FORTRAN_BLAS_STYLE_LOWCASE)

#define F77_sgetri sgetri
#define F77_dgetri dgetri
#define F77_sgetrf sgetrf
#define F77_dgetrf dgetrf
#define F77_strtrs strtrs
#define F77_dtrtrs dtrtrs
#define F77_stptrs stptrs
#define F77_dtptrs dtptrs
#define F77_slaswp slaswp
#define F77_dlaswp dlaswp

#else

// Nothing defined: should not happen

#error "Must define either LAMA_FOTRAN_BLAS_STYLE_UNDERSCORE, LAMA_FORTRAN_BLAS_STYLE_UPCASE, or LAMA_FORTRAN_BLAS_STYLE_LOWCASE"

#endif /* LAMA_FORTRAN_BLAS_STYLE */

void F77_sgetrf( const int* m, const int* n, float* a, const int* lda, int* ipivot, int* info );
void F77_dgetrf( const int* m, const int* n, double* a, const int* lda, int* ipivot, int* info );
void F77_sgetri( const int* n, float* a, const int* lda, int* ipivot, float* work, const int* ldwork, int* info );
void F77_dgetri( const int* n, double* a, const int* lda, int* ipivot, double* work, const int* ldwork, int* info );
void F77_stptrs(
    char* uplo,
    char* transa,
    char* diag,
    const int* n,
    const int* nrhs,
    const float* ap,
    float* b,
    const int* ldb,
    int* info );
void F77_dtptrs(
    char* uplo,
    char* transa,
    char* diag,
    const int* n,
    const int* nrhs,
    const double* ap,
    double* b,
    const int* ldb,
    int* info );
int F77_slaswp(
    const int* n,
    float* a,
    const int* lda,
    const int* k1,
    const int* k2,
    const int* ipiv,
    const int* incx );
int F77_dlaswp(
    const int* n,
    double* a,
    const int* lda,
    const int* k1,
    const int* k2,
    const int* ipiv,
    const int* incx );

#ifdef __cplusplus
} /*extern "C"*/
#endif /*__cplusplus*/


using common::unique_ptr;

namespace lama
{

/* ------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( LAPACK_LAPACK::logger, "LAPACK.LAPACK" )

/* ------------------------------------------------------------------------- */
/*      getrf<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACK_LAPACK::getrf(
    const CBLAS_ORDER order,
    const int m,
    const int n,
    float* const A,
    const int lda,
    int* const ipiv )
{
    LAMA_REGION( "LAPACK.LAPACK.getrf<float>" )

    LAMA_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

    int info = 0;

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    if( order == CblasColMajor )
    {
        F77_sgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
    }
    else if( m == n && n == lda )
    {
        for( int i = 0; i < m; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( A[i * n + j], A[j * m + i] );
            }
        }

        F77_sgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );

        for( int i = 0; i < m; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( A[i * n + j], A[j * m + i] );
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "row major only supported for square matrices" );
    }

    for( int i = 0; i < m; ++i )
    {
        --ipiv[i]; // Fortran numbering from 1 to n ->  0 to n-1
    }

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getrf<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
IndexType LAPACK_LAPACK::getrf(
    const CBLAS_ORDER order,
    const int m,
    const int n,
    double* const A,
    const int lda,
    int* const ipiv )
{
    LAMA_REGION( "LAPACK.LAPACK.getrf<double>" )

    LAMA_LOG_INFO( logger, "getrf<double> for A of size " << m << " x " << n )

    int info = 0;

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    if( order == CblasColMajor )
    {
        F77_dgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
    }
    else if( m == n && n == lda )
    {
        for( int i = 0; i < m; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( A[i * n + j], A[j * m + i] );
            }
        }

        F77_dgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );

        for( int i = 0; i < m; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( A[i * n + j], A[j * m + i] );
            }
        }
    }
    else
    {
        COMMON_THROWEXCEPTION( "row major only supported for square matrices" );
    }

    for( int i = 0; i < m; ++i )
    {
        --ipiv[i]; // Fortran numbering from 1 to n ->  0 to n-1
    }

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getinv<float>                                                        */
/* ------------------------------------------------------------------------- */

template<>
void LAPACK_LAPACK::getinv( const IndexType n, float* a, const IndexType lda )
{
    LAMA_REGION( "LAPACK.LAPACK.getinv<float>" )

    int info = 0;

    // unique_ptr, delete by destructor, also done in case of exception

    scoped_array<IndexType> ipiv( new IndexType[n] );

#ifdef F77_INT
    F77_INT F77_N = n, F77_lda = lda, F77_info = info;
#else
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    LAMA_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses Fortran interface" )

    F77_sgetrf( &F77_N, &F77_N, a, &F77_lda, ipiv.get(), &F77_info );

    if( F77_info )
    {
        COMMON_THROWEXCEPTION( "LAPACK sgetrf failed, info = " << F77_info )
    }

    scoped_array<float> work( new float[n] );

    F77_sgetri( &F77_N, a, &F77_lda, ipiv.get(), work.get(), &F77_N, &F77_info );

    if( F77_info )
    {
        COMMON_THROWEXCEPTION( "LAPACK sgetri failed, info = " << F77_info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getinv<double>                                                       */
/* ------------------------------------------------------------------------- */

template<>
void LAPACK_LAPACK::getinv( const IndexType n, double* a, const IndexType lda )
{
    LAMA_REGION( "LAPACK.LAPACK.getinv<double>" )

    int info = 0;

    scoped_array<IndexType> ipiv( new IndexType[n] );

#ifdef F77_INT
    F77_INT F77_N = n, F77_lda = lda, F77_info = info;
#else
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    LAMA_LOG_INFO( logger, "getinv<double> for " << n << " x " << n << " matrix, uses Fortran interface" )

    F77_dgetrf( &F77_N, &F77_N, a, &F77_lda, ipiv.get(), &F77_info );

    if( F77_info )
    {
        COMMON_THROWEXCEPTION( "LAPACK dgetrf failed, info = " << F77_info )
    }

    scoped_array<double> work( new double[n] );

    F77_dgetri( &F77_N, a, &F77_lda, ipiv.get(), work.get(), &F77_N, &F77_info );

    if( F77_info )
    {
        COMMON_THROWEXCEPTION( "LAPACK dgetri failed, info = " << F77_info )
    }
}

/* ------------------------------------------------------------------------- */
/*      getri<float>                                                         */
/* ------------------------------------------------------------------------- */

template<>
int LAPACK_LAPACK::getri( const CBLAS_ORDER order, const int n, float* const a, const int lda, int* const ipiv )
{
    LAMA_REGION( "LAPACK.LAPACK.getri<float>" )

    LAMA_LOG_INFO( logger, "getri<float> for A of size " << n << " x " << n )

    int info = 0;

    // translate C indexes into  Fortran Indexes for ipiv

    for( int i = 0; i < n; ++i )
    {
        ++ipiv[i];
    }

    // transpose if not column major order

    if( order != CblasColMajor )
    {
        LAMA_ASSERT_EQUAL_ERROR( lda, n )

        for( int i = 0; i < n; ++i )
        {
            // swap row and column

            for( int j = i + 1; j < n; ++j )
            {
                std::swap( a[i * n + j], a[j * n + i] );
            }
        }
    }

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    scoped_array<float> work( new float[n] );

    F77_sgetri( &F77_N, a, &F77_lda, ipiv, work.get(), &F77_N, &F77_info );

    if( order != CblasColMajor )
    {
        // transpose back

        for( int i = 0; i < n; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( a[i * n + j], a[j * n + i] );
            }
        }
    }

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* ------------------------------------------------------------------------- */
/*      getri<double>                                                        */
/* ------------------------------------------------------------------------- */

template<>
int LAPACK_LAPACK::getri( const CBLAS_ORDER order, const int n, double* const a, const int lda, int* const ipiv )
{
    LAMA_REGION( "LAPACK.LAPACK.getri<double>" )

    LAMA_LOG_INFO( logger, "getri<double> for A of size " << n << " x " << n )

    int info = 0;

    // translate C indexes into  Fortran Indexes for ipiv

    for( int i = 0; i < n; ++i )
    {
        ++ipiv[i];
    }

    // transpose if not column major order

    if( order != CblasColMajor )
    {
        LAMA_ASSERT_EQUAL_ERROR( lda, n )

        for( int i = 0; i < n; ++i )
        {
            // swap row and column

            for( int j = i + 1; j < n; ++j )
            {
                std::swap( a[i * n + j], a[j * n + i] );
            }
        }
    }

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    scoped_array<double> work( new double[n] );

    F77_dgetri( &F77_N, a, &F77_lda, ipiv, work.get(), &F77_N, &F77_info );

    if( order != CblasColMajor )
    {
        // transpose back

        for( int i = 0; i < n; ++i )
        {
            for( int j = i + 1; j < n; ++j )
            {
                std::swap( a[i * n + j], a[j * n + i] );
            }
        }
    }

    if( info < 0 )
    {
        COMMON_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if( info > 0 )
    {
        COMMON_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<>
int LAPACK_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const float* AP,
    float* B,
    const int ldb )
{
    LAMA_REGION( "LAPACK.LAPACK.tptrs<float>" )

    int info = 0;

    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

#ifdef F77_CHAR
    F77_CHAR F77_UL, F77_TA, F77_DI;
#else
#   define F77_UL &UL
#   define F77_TA &TA
#   define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_n = n, F77_nrhs = nrhs, F77_ldb = ldb;
#else
#   define F77_n n
#   define F77_nrhs nrhs
#   define F77_ldb ldb
#endif

#ifdef F77_CHAR
    F77_UL = C2F_CHAR( &UL );
    F77_TA = C2F_CHAR( &TA );
    F77_DI = C2F_CHAR( &DI );
#endif

    LAMA_LOG_INFO( logger,
                   "tptrs<float>, n = " << n << ", nrhs = " << nrhs << ", order = " << order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    if( order == CblasColMajor )
    {
        F77_stptrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, AP, B, &F77_ldb, &info );
    }
    else if( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        COMMON_THROWEXCEPTION( "row major order not supported for tptrs" );
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<>
int LAPACK_LAPACK::tptrs(
    const CBLAS_ORDER order,
    const CBLAS_UPLO uplo,
    const CBLAS_TRANSPOSE trans,
    const CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const double* AP,
    double* B,
    const int ldb )
{
    LAMA_REGION( "LAPACK.LAPACK.tptrs<double>" )

    int info = 0;

    char UL = BLASHelper::lapack_uplo( uplo );
    char TA = BLASHelper::lapack_transpose( trans );
    char DI = BLASHelper::lapack_diag( diag );

#ifdef F77_CHAR
    F77_CHAR F77_UL, F77_TA, F77_DI;
#else
#   define F77_UL &UL
#   define F77_TA &TA
#   define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_n = n, F77_nrhs = nrhs, F77_ldb = ldb;
#else
#   define F77_n n
#   define F77_nrhs nrhs
#   define F77_ldb ldb
#endif

#ifdef F77_CHAR
    F77_UL = C2F_CHAR( &UL );
    F77_TA = C2F_CHAR( &TA );
    F77_DI = C2F_CHAR( &DI );
#endif

    LAMA_LOG_INFO( logger,
                   "tptrs<double>, n = " << n << ", nrhs = " << nrhs << ", order = " << order << ", UL = " << UL << ", TA = " << TA << ", DI = " << DI );

    if( order == CblasColMajor )
    {
        F77_dtptrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, AP, B, &F77_ldb, &info );
    }
    else if( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        COMMON_THROWEXCEPTION( "row major order not supported for tptrs" );
    }

    return info;
}

/* --------------------------------------------------------------------------- */

template<>
void LAPACK_LAPACK::laswp(
    const CBLAS_ORDER order,
    const int N,
    float* A,
    const int LDA,
    const int K1,
    const int K2,
    const int* ipiv,
    const int INCX,
    SyncToken* syncToken )
{
    LAMA_REGION( "LAPACK.LAPACK.laswp<float>" )

    int i = K1;

    if( order == CblasRowMajor )
    {
        for( i = K1; i < K2; ++i )
        {
            if( ipiv[i * INCX] == i )
            {
                continue;
            }

            BLAS_BLAS1::swap<float>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX, syncToken );
        }
    }
    else if( order == CblasColMajor )
    {
#ifdef F77_INT
        F77_INT F77_N = N, F77_LDA = LDA, F77_K1 = K1, F77_K2 = K2, F77_INCX = INCX;
        F77_INT* F77_IPIV = ipiv;
#else
#   ifdef F77_N
#       undef F77_N
#   endif
#   define F77_N    N
#   define F77_LDA  LDA
#   define F77_K1   K1
#   define F77_K2   K2
#   define F77_INCX INCX
#   define F77_IPIV ipiv
#endif
        F77_slaswp( &F77_N, A, &F77_LDA, &F77_K1, &F77_K2, F77_IPIV, &F77_INCX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( 0, 1, "cblas_slaswp", "Illegal order setting, %d\n", order );
    }
}

/* --------------------------------------------------------------------------- */

template<>
void LAPACK_LAPACK::laswp(
    const CBLAS_ORDER order,
    const int N,
    double* A,
    const int LDA,
    const int K1,
    const int K2,
    const int* ipiv,
    const int INCX,
    SyncToken* syncToken )
{
    LAMA_REGION( "LAPACK.LAPACK.laswp<double>" )

    int i = K1;

    if( order == CblasRowMajor )
    {
        for( i = K1; i < K2; ++i )
        {
            if( ipiv[i * INCX] == i )
            {
                continue;
            }

            BLAS_BLAS1::swap<double>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX, syncToken );
        }
    }
    else if( order == CblasColMajor )
    {
#ifdef F77_INT
        F77_INT F77_N = N, F77_LDA = LDA, F77_K1 = K1, F77_K2 = K2, F77_INCX = INCX;
        F77_INT* F77_IPIV = ipiv;
#else
#   define F77_N    N
#   define F77_LDA  LDA
#   define F77_K1   K1
#   define F77_K2   K2
#   define F77_INCX INCX
#   define F77_IPIV ipiv
#endif
        F77_dlaswp( &F77_N, A, &F77_LDA, &F77_K1, &F77_K2, F77_IPIV, &F77_INCX );
    }
    else
    {
        COMMON_THROWEXCEPTION( "Illegal order setting: " << order );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void LAPACK_LAPACK::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions
    //       ( e.g. routine CUDABLAS1::sum<ValueType> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, getrf, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getrf, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, getri, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getri, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, getinv, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, getinv, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, tptrs, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, tptrs, double )

    LAMA_INTERFACE_REGISTER_T( BLAS, laswp, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, laswp, double )

    // other routines are not used by LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the LAPACK routines                               */
/* --------------------------------------------------------------------------- */

bool LAPACK_LAPACK::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool LAPACK_LAPACK::initialized = registerInterface();

} /* namespace lama */
