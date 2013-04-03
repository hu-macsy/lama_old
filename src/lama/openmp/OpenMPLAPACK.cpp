/**
 * @file OpenMPLAPACK.cpp
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
 * @brief OpenMPLAPACK.cpp
 * @author lschubert
 * @date 02.07.2012
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPLAPACK.hpp>

// others
#include <lama/openmp/OpenMPBLAS1.hpp>

#include <lama/BLASInterface.hpp>
#include <boost/scoped_array.hpp>

// macros
#include <lama/macros/unused.hpp>

#ifdef __cplusplus
extern "C"
{
#endif /*__cplusplus*/

#ifdef LAMA_HAVE_MKL

#   include <mkl_lapack.h>
#   include <mkl_lapacke.h>

#define F77_sgetrf sgetrf_
#define F77_dgetrf dgetrf_
#define F77_sgetri sgetri_
#define F77_dgetri dgetri_
#define F77_strtrs strtrs_
#define F77_dtrtrs dtrtrs_
#define F77_stptrs stptrs_
#define F77_dtptrs dtptrs_

#else /* LAMA_BLAS_MKL */

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
    void F77_strtrs(
        char* uplo,
        char* transa,
        char* diag,
        const int* n,
        const int* nrhs,
        const float* a,
        const int* lda,
        float* b,
        const int* ldb,
        int* info );
    void F77_dtrtrs(
        char* uplo,
        char* transa,
        char* diag,
        const int* n,
        const int* nrhs,
        const double* a,
        const int* lda,
        double* b,
        const int* ldb,
        int* info );
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

#endif /* LAMA_BLAS_MKL */

#ifdef __cplusplus
} /*extern "C"*/
#endif /*__cplusplus*/

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPLAPACK::logger, "OpenMP.LAPACK" )

template<>
IndexType OpenMPLAPACK::getrf(
    const enum CBLAS_ORDER order,
    const int m,
    const int n,
    float* const A,
    const int lda,
    int* const ipiv )
{
    LAMA_LOG_INFO( logger, "getrf<float> for A of size " << m << " x " << n )

    int info = 0;

#ifdef LAMA_HAVE_MKL

    LAMA_LOG_INFO( logger, "will use MKL routine" )

    if ( order == CblasColMajor )
    {
        info = LAPACKE_sgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );
    }
    else
    {
//        printf( "Calling LAPACKE_sgetrf( %d, %d, %d, %f, %d, %d )\n,", LAPACK_COL_MAJOR, m, n, A[0], lda, ipiv[0] );
        info = LAPACKE_sgetrf( LAPACK_ROW_MAJOR, m, n, A, lda, ipiv );
    }

#else //LAMA_HAVE_MKL
    int i = 0;
#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    if ( order == CblasColMajor )
    {
        F77_sgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
    }
    else
    {
        if ( m == n && n == lda )
        {
            float tmp = 0.0f;

            for ( i = 0; i < m; ++i )
            {
                int j = 0;

                // swap row and column
                for ( j = i + 1; j < n; ++j )
                {
                    tmp = A[i * n + j];
                    A[i * n + j] = A[j * m + i];
                    A[j + m + i] = tmp;
                }
            }
        }
        else
        {
            //TODO: Transpose non-square or blocked intput matrix
            return -1 * ( m + n + lda );
        }

        F77_sgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
    }

    //Fix Fortran Indexes in ipiv array
    #pragma omp parallel for

    for ( i = 0; i < m; ++i )
    {
        --ipiv[i];
    }

#endif //LAMA_HAVE_MKL
    if ( info < 0 )
    {
        LAMA_THROWEXCEPTION( "illegal argument " << ( -info ) )
    }
    else if ( info > 0 )
    {
        LAMA_THROWEXCEPTION( "value(" << info << "," << info << ")" << " is exactly zero" )
    }

    return info;
}

template<>
IndexType OpenMPLAPACK::getrf(
    const enum CBLAS_ORDER order,
    const int m,
    const int n,
    double* const A,
    const int lda,
    int* const ipiv )
{
    int info = 0;
#ifdef LAMA_HAVE_MKL

    if ( order == CblasColMajor )
    {
        info = LAPACKE_dgetrf( LAPACK_COL_MAJOR, m, n, A, lda, ipiv );
    }
    else
    {
        info = LAPACKE_dgetrf( LAPACK_ROW_MAJOR, m, n, A, lda, ipiv );
    }

#else //LAMA_HAVE_MKL
    int i = 0;
#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_info = info;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    if ( order == CblasColMajor )
    {
        F77_dgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
    }
    else
    {
        //TODO: Transpose intput matrix
        return -1 * ( m + n + lda );
        //F77_dgetrf( &F77_M, &F77_N, A, &F77_lda, ipiv, &F77_info );
    }

    #pragma omp parallel for

    //Fix Fortran Indexes in ipiv array
    for ( i = 0; i < m; ++i )
    {
        --ipiv[i];
    }

#endif //LAMA_HAVE_MKL
    return info;
}

template<>
void OpenMPLAPACK::getinv( const IndexType n, float* a, const IndexType lda )
{
    int info = 0;

    // scoped array, will also be freed in case of exception

    boost::scoped_array<IndexType> ipiv( new IndexType[n] );

#ifdef LAMA_HAVE_MKL

    LAMA_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKE_sgetrf( LAPACK_COL_MAJOR, n, n, a, lda, ipiv.get() );

    // return error if factorization did not work

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL sgetrf failed, info = " << info )
    }

    info = LAPACKE_sgetri( LAPACK_COL_MAJOR, n, a, lda, ipiv.get() );

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL sgetri failed, info = " << info )
    }

#else //LAMA_HAVE_MKL
#ifdef F77_INT
    F77_INT F77_N = n, F77_lda = lda, F77_info = info;
#else
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    LAMA_LOG_INFO( logger, "getinv<float> for " << n << " x " << n << " matrix, uses Fortran interface" )

    F77_sgetrf( &F77_N, &F77_N, a, &F77_lda, ipiv.get(), &F77_info );

    if ( F77_info )
    {
        LAMA_THROWEXCEPTION( "LAPACK sgetrf failed, info = " << F77_info )
    }

    boost::scoped_array<float> work( new float[n] );

    F77_sgetri( &F77_N, a, &F77_lda, ipiv.get(), work.get(), &F77_N, &F77_info );

    if ( F77_info )
    {
        LAMA_THROWEXCEPTION( "LAPACK sgetri failed, info = " << F77_info )
    }

#endif // LAMA_HAVE_MKL
}

template<>
void OpenMPLAPACK::getinv( const IndexType n, double* a, const IndexType lda )
{
    int info = 0;

    boost::scoped_array<IndexType> ipiv( new IndexType[n] );

#ifdef LAMA_HAVE_MKL

    LAMA_LOG_INFO( logger, "getinv<double> for " << n << " x " << n << " matrix, uses MKL" )

    info = LAPACKE_dgetrf( LAPACK_COL_MAJOR, n, n, a, lda, ipiv.get() );

    // return error if factorization did not work

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL dgetrf failed, info = " << info )
    }

    info = LAPACKE_dgetri( LAPACK_COL_MAJOR, n, a, lda, ipiv.get() );

    if ( info )
    {
        LAMA_THROWEXCEPTION( "MKL dgetri failed, info = " << info )
    }

#else //LAMA_HAVE_MKL
#ifdef F77_INT
    F77_INT F77_N = n, F77_lda = lda, F77_info = info;
#else
#define F77_N n
#define F77_lda lda
#define F77_info info
#endif

    LAMA_LOG_INFO( logger, "getinv<double> for " << n << " x " << n << " matrix, uses Fortran interface" )

    F77_dgetrf( &F77_N, &F77_N, a, &F77_lda, ipiv.get(), &F77_info );

    if ( F77_info )
    {
        LAMA_THROWEXCEPTION( "LAPACK dgetrf failed, info = " << F77_info )
    }

    boost::scoped_array<double> work( new double[n] );

    F77_dgetri( &F77_N, a, &F77_lda, ipiv.get(), work.get(), &F77_N, &F77_info );

    if ( F77_info )
    {
        LAMA_THROWEXCEPTION( "LAPACK dgetri failed, info = " << F77_info )
    }

#endif // LAMA_HAVE_MKL
}

#ifdef LAMA_HAVE_MKL
template<>
int OpenMPLAPACK::getri(
    const enum CBLAS_ORDER order,
    const int n,
    float* const a,
    const int lda,
    int* const ipiv )
{
    int info = 0;
    info = LAPACKE_sgetri( order, n, a, lda, ipiv );
    return info;
}
#else
template<>
int OpenMPLAPACK::getri(
    const enum CBLAS_ORDER UNUSED( order ),
    const int UNUSED( n ),
    float* const UNUSED( a ),
    const int UNUSED( lda ),
    int* const UNUSED( ipiv ) )
{
    LAMA_THROWEXCEPTION( "Not implemented." )
    return -1;
}
#endif

#ifdef LAMA_HAVE_MKL
template<>
int OpenMPLAPACK::getri(
    const enum CBLAS_ORDER order,
    const int n,
    double* const a,
    const int lda,
    int* const ipiv )
{
    int info = 0;
    info = LAPACKE_dgetri( order, n, a, lda, ipiv );
    return info;
}
#else
template<>
int OpenMPLAPACK::getri(
    const enum CBLAS_ORDER UNUSED( order ),
    const int UNUSED( n ),
    double* const UNUSED( a ),
    const int UNUSED( lda ),
    int* const UNUSED( ipiv ) )
{
    LAMA_THROWEXCEPTION( "Not implemented." )
    return -1;
}
#endif

template<>
int OpenMPLAPACK::trtrs(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const float* A,
    const int lda,
    float* B,
    const int ldb )
{
    char* f_name = (char*) "strtrs";
    int info = 0;
    char UL, TA, DI;
#ifdef F77_CHAR
    F77_CHAR F77_UL, F77_TA, F77_DI;
#else
#   define F77_UL &UL
#   define F77_TA &TA
#   define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_n = n, F77_nrhs = nrhs, F77_lda = lda, F77_ldb = ldb;
#else
#   define F77_n n
#   define F77_nrhs nrhs
#   define F77_lda lda
#   define F77_ldb ldb
#endif

    switch ( uplo )
    {
    case CblasUpper:
        UL = 'U';
        break;
    case CblasLower:
        UL = 'L';
        break;
    default:
        info = 2;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal uplo setting." );
        UL = 'U';
    }

    switch ( trans )
    {
    case CblasNoTrans:
        TA = 'N';
        break;
    case CblasTrans:
        TA = 'T';
        break;
    case CblasConjTrans:
        TA = 'C';
        break;
    default:
        info = 3;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal trans setting." );
        TA = 'N';
    }

    switch ( diag )
    {
    case CblasNonUnit:
        DI = 'N';
        break;
    case CblasUnit:
        DI = 'U';
        break;
    default:
        info = 4;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal diag setting." );
        DI = 'N';
    }

#ifdef F77_CHAR
    F77_UL = C2F_CHAR( &UL );
    F77_TA = C2F_CHAR( &TA );
    F77_DI = C2F_CHAR( &DI );
#endif
#ifdef LAMA_HAVE_MKL

    if ( order == CblasColMajor )
    {
        info = LAPACKE_strtrs( LAPACK_COL_MAJOR, UL, TA, DI, n, nrhs, A, lda, B, ldb );
    }
    else if ( order == CblasRowMajor )
    {
        info = LAPACKE_strtrs( LAPACK_ROW_MAJOR, UL, TA, DI, n, nrhs, A, lda, B, ldb );
    }
    else
    {
        info = 1;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal order setting." );
    }

#else /* LAMA_HAVE_MKL */

    if ( order == CblasColMajor )
    {
        F77_strtrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, A, &F77_lda, B, &F77_ldb, &info );
    }
    else if ( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        info = -1 * ( n + nrhs + lda + ldb );
    }

#endif/* LAMA_HAVE_MKL */
    return info;
}

template<>
int OpenMPLAPACK::trtrs(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const double* A,
    const int lda,
    double* B,
    const int ldb )
{
    char* f_name = (char*) "dtrtrs";
    int info = 0;
    char UL, TA, DI;
#ifdef F77_CHAR
    F77_CHAR F77_UL, F77_TA, F77_DI;
#else
#   define F77_UL &UL
#   define F77_TA &TA
#   define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_n = n, F77_nrhs = nrhs, F77_lda = lda, F77_ldb = ldb;
#else
#   define F77_n n
#   define F77_nrhs nrhs
#   define F77_lda lda
#   define F77_ldb ldb
#endif

    switch ( uplo )
    {
    case CblasUpper:
        UL = 'U';
        break;
    case CblasLower:
        UL = 'L';
        break;
    default:
        info = 2;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal uplo setting." );
        UL = 'U';
    }

    switch ( trans )
    {
    case CblasNoTrans:
        TA = 'N';
        break;
    case CblasTrans:
        TA = 'T';
        break;
    case CblasConjTrans:
        TA = 'C';
        break;
    default:
        info = 3;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal trans setting." );
        TA = 'N';
    }

    switch ( diag )
    {
    case CblasNonUnit:
        DI = 'N';
        break;
    case CblasUnit:
        DI = 'U';
        break;
    default:
        info = 4;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal diag setting." );
        DI = 'N';
    }

#ifdef F77_CHAR
    F77_UL = C2F_CHAR( &UL );
    F77_TA = C2F_CHAR( &TA );
    F77_DI = C2F_CHAR( &DI );
#endif
#ifdef LAMA_HAVE_MKL

    if ( order == CblasColMajor )
    {
        info = LAPACKE_dtrtrs( LAPACK_COL_MAJOR, UL, TA, DI, n, nrhs, A, lda, B, ldb );
    }
    else if ( order == CblasRowMajor )
    {
        info = LAPACKE_dtrtrs( LAPACK_ROW_MAJOR, UL, TA, DI, n, nrhs, A, lda, B, ldb );
    }
    else
    {
        info = 1;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal order setting." );
    }

#else

    if ( order == CblasColMajor )
    {
        F77_dtrtrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, A, &F77_lda, B, &F77_ldb, &info );
    }
    else if ( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        info = -1 * ( n + nrhs + lda + ldb );
    }

#endif
    return info;
}

template<>
int OpenMPLAPACK::tptrs(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const float* AP,
    float* B,
    const int ldb )
{
    char* f_name = (char*) "stptrs";
    int info = 0;
    char UL, TA, DI;
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

    switch ( uplo )
    {
    case CblasUpper:
        UL = 'U';
        break;
    case CblasLower:
        UL = 'L';
        break;
    default:
        info = 2;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal uplo setting." );
        UL = 'U';
    }

    switch ( trans )
    {
    case CblasNoTrans:
        TA = 'N';
        break;
    case CblasTrans:
        TA = 'T';
        break;
    case CblasConjTrans:
        TA = 'C';
        break;
    default:
        info = 3;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal trans setting." );
        TA = 'N';
    }

    switch ( diag )
    {
    case CblasNonUnit:
        DI = 'N';
        break;
    case CblasUnit:
        DI = 'U';
        break;
    default:
        info = 4;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal diag setting." );
        DI = 'N';
    }

#ifdef F77_CHAR
    F77_UL = C2F_CHAR( &UL );
    F77_TA = C2F_CHAR( &TA );
    F77_DI = C2F_CHAR( &DI );
#endif
#ifdef LAMA_HAVE_MKL

    if ( order == CblasColMajor )
    {
        info = LAPACKE_stptrs( LAPACK_COL_MAJOR, UL, TA, DI, n, nrhs, AP, B, ldb );
    }
    else if ( order == CblasRowMajor )
    {
        info = LAPACKE_stptrs( LAPACK_ROW_MAJOR, UL, TA, DI, n, nrhs, AP, B, ldb );
    }
    else
    {
        info = 1;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal order setting." );
    }

#else

    if ( order == CblasColMajor )
    {
        F77_stptrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, AP, B, &F77_ldb, &info );
    }
    else if ( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        info = -1 * ( n + nrhs + ldb );
    }

#endif
    return info;
}

template<>
int OpenMPLAPACK::tptrs(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO uplo,
    const enum CBLAS_TRANSPOSE trans,
    const enum CBLAS_DIAG diag,
    const int n,
    const int nrhs,
    const double* AP,
    double* B,
    const int ldb )
{
    char* f_name = (char*) "dtptrs";
    int info = 0;
    char UL, TA, DI;
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

    switch ( uplo )
    {
    case CblasUpper:
        UL = 'U';
        break;
    case CblasLower:
        UL = 'L';
        break;
    default:
        info = 2;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal uplo setting." );
        UL = 'U';
    }

    switch ( trans )
    {
    case CblasNoTrans:
        TA = 'N';
        break;
    case CblasTrans:
        TA = 'T';
        break;
    case CblasConjTrans:
        TA = 'C';
        break;
    default:
        info = 3;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal trans setting." );
        TA = 'N';
    }

    switch ( diag )
    {
    case CblasNonUnit:
        DI = 'N';
        break;
    case CblasUnit:
        DI = 'U';
        break;
    default:
        info = 4;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal diag setting." );
        DI = 'N';
    }

#ifdef F77_CHAR
    F77_UL = C2F_CHAR( &UL );
    F77_TA = C2F_CHAR( &TA );
    F77_DI = C2F_CHAR( &DI );
#endif
#ifdef LAMA_HAVE_MKL

    if ( order == CblasColMajor )
    {
        info = LAPACKE_dtptrs( LAPACK_COL_MAJOR, UL, TA, DI, n, nrhs, AP, B, ldb );
    }
    else if ( order == CblasRowMajor )
    {
        info = LAPACKE_dtptrs( LAPACK_ROW_MAJOR, UL, TA, DI, n, nrhs, AP, B, ldb );
    }
    else
    {
        info = 1;
        BLASHelper::XERBLA_cpu( 0, info, f_name, "Illegal order setting." );
    }

#else

    if ( order == CblasColMajor )
    {
        F77_dtptrs( F77_UL, F77_TA, F77_DI, &F77_n, &F77_nrhs, AP, B, &F77_ldb, &info );
    }
    else if ( order == CblasRowMajor )
    {
        // TODO: transpose matrix.
        info = -1 * ( n + nrhs + ldb );
    }

#endif
    return info;
}

template<>
void OpenMPLAPACK::laswp(
    const enum CBLAS_ORDER order,
    const int N,
    float* A,
    const int LDA,
    const int K1,
    const int K2,
    const int* ipiv,
    const int INCX,
    SyncToken* syncToken )
{
    int i = K1;

    if ( order == CblasRowMajor )
    {
        for ( i = K1; i < K2; ++i )
        {
            if ( ipiv[i * INCX] == i )
            {
                continue;
            }

            OpenMPBLAS1::swap<float>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX, syncToken );
        }
    }
    else if ( order == CblasColMajor )
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

template<>
void OpenMPLAPACK::laswp(
    const enum CBLAS_ORDER order,
    const int N,
    double* A,
    const int LDA,
    const int K1,
    const int K2,
    const int* ipiv,
    const int INCX,
    SyncToken* syncToken )
{
    int i = K1;

    if ( order == CblasRowMajor )
    {
        for ( i = K1; i < K2; ++i )
        {
            if ( ipiv[i * INCX] == i )
            {
                continue;
            }

            OpenMPBLAS1::swap<double>( N, &A[ipiv[i * INCX] * LDA], INCX, &A[i * LDA], INCX, syncToken );
        }
    }
    else if ( order == CblasColMajor )
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
        BLASHelper::XERBLA_cpu( 0, 1, "cblas_dlaswp", "Illegal order setting, %d\n", order );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPLAPACK::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

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

} /* namespace lama */
