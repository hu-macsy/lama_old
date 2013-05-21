/**
 * @file OpenMPBLAS2.cpp
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
 * @brief OpenMPBLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPBLAS2.hpp>
#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPBLAS2::logger, "OpenMP.BLAS2" )

/** gemv */

template<>
void OpenMPBLAS2::gemv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const IndexType M,
    const IndexType N,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* X,
    const IndexType incX,
    const float beta,
    float* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "gemv<float>: M = " << M << ", N = " << N << ", LDA = " << lda 
                    << ", incX = " << incX << ", incY = " << incY 
                    << ", alpha = " << alpha << ", beta = " << beta )

    if ( N == 0 )
    {
        return;   // empty X, Y, A  
    }

    // M == 0: empty A, we can handle this here

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
#ifdef F77_CHAR
    F77_CHAR F77_TA;
#else
#define F77_TA &TA
#endif
#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_M M
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_sgemv( F77_TA, &F77_M, &F77_N, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_sgemv( F77_TA, &F77_N, &F77_M, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sgemv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::gemv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const IndexType M,
    const IndexType N,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* X,
    const IndexType incX,
    const double beta,
    double* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "gemv<double>: M = " << M << ", N = " << N << ", LDA = " << lda 
                    << ", incX = " << incX << ", incY = " << incY 
                    << ", alpha = " << alpha << ", beta = " << beta )

    if ( N == 0 )
    {
        return;   // empty X, Y, A  
    }

    // M == 0: empty A, we can handle this here

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
#ifdef F77_CHAR
    F77_CHAR F77_TA;
#else
#define F77_TA &TA
#endif
#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_M M
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgemv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_dgemv( F77_TA, &F77_M, &F77_N, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgemv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_dgemv( F77_TA, &F77_N, &F77_M, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dgemv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** symv */

template<>
void OpenMPBLAS2::symv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* X,
    const IndexType incX,
    const float beta,
    float* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "symv<float>: N = " << N << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssymv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssymv( F77_UL, &F77_N, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssymv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssymv( F77_UL, &F77_N, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssymv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
    // ssymv( &uplo, &n, &alpha, A, &lda, x, &incX, &beta, y, &incY);
}

template<>
void OpenMPBLAS2::symv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* X,
    const IndexType incX,
    const double beta,
    double* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "symv<double>: N = " << N << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsymv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsymv( F77_UL, &F77_N, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsymv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsymv( F77_UL, &F77_N, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsymv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** trmv */

template<>
void OpenMPBLAS2::trmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const float* A,
    const IndexType lda,
    float* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "symv<float>: N = " << N << ", LDA = " << lda << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_incX = incX;
#else
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strmv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strmv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_strmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::trmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const double* A,
    const IndexType lda,
    double* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "symv<double>: N = " << N << ", LDA = " << lda << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_incX = incX;
#else
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrmv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrmv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtrmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** trsv */

template<>
void OpenMPBLAS2::trsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const float* A,
    const IndexType lda,
    float* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "trsv<float>: N = " << N << ", LDA = " << lda << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_incX = incX;
#else
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strsv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strsv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_strsv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::trsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const double* A,
    const IndexType lda,
    double* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "trsv<doubel>: N = " << N << ", LDA = " << lda << ", incX = " << incX )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_incX = incX;
#else
#define F77_N N
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrsv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrsv( F77_UL, F77_TA, F77_DI, &F77_N, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtrsv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** gbmv */

template<>
void OpenMPBLAS2::gbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const IndexType M,
    const IndexType N,
    const IndexType KL,
    const IndexType KU,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* X,
    const IndexType incX,
    const float beta,
    float* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "gbmv<float>: M = " << M << ", N = " << N << ", KL = " << KL << ", KU = " << KU << ", alpha = " << alpha << ", LDA = " << lda << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
#ifdef F77_CHAR
    F77_CHAR F77_TA;
#else
#define F77_TA &TA
#endif
#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
    F77_INT F77_KL = KL, F77_KU = KU;
#else
#define F77_M M
#define F77_N N
#define F77_lda lda
#define F77_KL KL
#define F77_KU KU
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_sgbmv( F77_TA, &F77_M, &F77_N, &F77_KL, &F77_KU, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_sgbmv( F77_TA, &F77_N, &F77_M, &F77_KU, &F77_KL, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sgbmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::gbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const IndexType M,
    const IndexType N,
    const IndexType KL,
    const IndexType KU,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* X,
    const IndexType incX,
    const double beta,
    double* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "gbmv<double>: M = " << M << ", N = " << N << ", KL = " << KL << ", KU = " << KU << ", alpha = " << alpha << ", LDA = " << lda << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
#ifdef F77_CHAR
    F77_CHAR F77_TA;
#else
#define F77_TA &TA
#endif
#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
    F77_INT F77_KL = KL, F77_KU = KU;
#else
#define F77_M M
#define F77_N N
#define F77_lda lda
#define F77_KL KL
#define F77_KU KU
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_dgbmv( F77_TA, &F77_M, &F77_N, &F77_KL, &F77_KU, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
#endif
        F77_dgbmv( F77_TA, &F77_N, &F77_M, &F77_KU, &F77_KL, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dgbmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
}

/** sbmv */

template<>
void OpenMPBLAS2::sbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const IndexType K,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* X,
    const IndexType incX,
    const float beta,
    float* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "sbmv<float>: N = " << N << ", K = " << K << ", alpha = " << alpha << ", LDA = " << lda << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_K = K, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssbmv( F77_UL, &F77_N, &F77_K, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssbmv( F77_UL, &F77_N, &F77_K, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssbmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::sbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const IndexType K,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* X,
    const IndexType incX,
    const double beta,
    double* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "sbmv<double>: N = " << N << ", K = " << K << ", alpha = " << alpha << ", LDA = " << lda << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_K = K, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsbmv( F77_UL, &F77_N, &F77_K, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsbmv( F77_UL, &F77_N, &F77_K, &alpha, A, &F77_lda, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsbmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** tbmv */

template<>
void OpenMPBLAS2::tbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const IndexType K,
    const float* A,
    const IndexType lda,
    float* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tbmv<float>: N = " << N << ", K = " << K << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_K = K, F77_incX = incX;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stbmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stbmv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stbmv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_stbmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::tbmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const IndexType K,
    const double* A,
    const IndexType lda,
    double* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tbmv<double>: N = " << N << ", K = " << K << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_K = K, F77_incX = incX;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtbmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtbmv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtbmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtbmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtbmv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtbmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
}

/** tbsv */

template<>
void OpenMPBLAS2::tbsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const IndexType K,
    const float* A,
    const IndexType lda,
    float* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tbsv<float>: N = " << N << ", K = " << K << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_K = K, F77_incX = incX;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stbsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stbsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stbsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stbsv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stbsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stbsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stbsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stbsv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_stbsv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::tbsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const IndexType K,
    const double* A,
    const IndexType lda,
    double* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tbsv<double>: N = " << N << ", K = " << K << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_lda = lda, F77_K = K, F77_incX = incX;
#else
#define F77_N N
#define F77_K K
#define F77_lda lda
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtbsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtbsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtbsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtbsv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtbsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtbsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtbsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtbsv( F77_UL, F77_TA, F77_DI, &F77_N, &F77_K, A, &F77_lda, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtbsv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** ger */

template<>
void OpenMPBLAS2::ger(
    const enum CBLAS_ORDER order,
    const IndexType M,
    const IndexType N,
    const float alpha,
    const float* X,
    const IndexType incX,
    const float* Y,
    const IndexType incY,
    float* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "ger<float>: M = " << M << ", N = " << N << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_M M
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#define F77_lda lda
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        F77_sger( &F77_M, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;
        F77_sger( &F77_N, &F77_M, &alpha, Y, &F77_incY, X, &F77_incX, A, &F77_lda );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sger", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::ger(
    const enum CBLAS_ORDER order,
    const IndexType M,
    const IndexType N,
    const double alpha,
    const double* X,
    const IndexType incX,
    const double* Y,
    const IndexType incY,
    double* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "ger<double>: M = " << M << ", N = " << N << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

#ifdef F77_INT
    F77_INT F77_M = M, F77_N = N, F77_lda = lda, F77_incX = incX, F77_incY = incY;
#else
#define F77_M M
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#define F77_lda lda
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        F77_dger( &F77_M, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;
        F77_dger( &F77_N, &F77_M, &alpha, Y, &F77_incY, X, &F77_incX, A, &F77_lda );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dger", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** syr */

template<>
void OpenMPBLAS2::syr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const float alpha,
    const float* X,
    const IndexType incX,
    float* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "syr<float>: N = " << N << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_lda = lda;
#else
#define F77_N N
#define F77_incX incX
#define F77_lda  lda
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssyr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssyr( F77_UL, &F77_N, &alpha, X, &F77_incX, A, &F77_lda );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssyr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssyr( F77_UL, &F77_N, &alpha, X, &F77_incX, A, &F77_lda );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssyr", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::syr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const double alpha,
    const double* X,
    const IndexType incX,
    double* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "syr<double>: N = " << N << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_lda = lda;
#else
#define F77_N N
#define F77_incX incX
#define F77_lda  lda
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsyr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsyr( F77_UL, &F77_N, &alpha, X, &F77_incX, A, &F77_lda );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsyr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsyr( F77_UL, &F77_N, &alpha, X, &F77_incX, A, &F77_lda );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsyr", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** syr2 */

template<>
void OpenMPBLAS2::syr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const float alpha,
    const float* X,
    const IndexType incX,
    const float* Y,
    const IndexType incY,
    float* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "syr2<float>: N = " << N << ", alpha = " << alpha << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_incY = incY, F77__lda = lda;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#define F77_lda  lda
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssyr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssyr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssyr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_ssyr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssyr2", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::syr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const double alpha,
    const double* X,
    const IndexType incX,
    const double* Y,
    const IndexType incY,
    double* A,
    const IndexType lda,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "syr2<double>: N = " << N << ", alpha = " << alpha << ", LDA = " << lda )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_incY = incY, F77__lda = lda;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#define F77_lda  lda
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsyr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsyr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsyr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dsyr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A, &F77_lda );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsyr2", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** spmv */

template<>
void OpenMPBLAS2::spmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const float alpha,
    const float* AP,
    const float* X,
    const IndexType incX,
    const float beta,
    float* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "spmv<float>: N = " << N << ", alpha = " << alpha << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sspmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_sspmv( F77_UL, &F77_N, &alpha, AP, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sspmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_sspmv( F77_UL, &F77_N, &alpha, AP, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sspmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
}

template<>
void OpenMPBLAS2::spmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const double alpha,
    const double* AP,
    const double* X,
    const IndexType incX,
    const double beta,
    double* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "spmv<double>: N = " << N << ", alpha = " << alpha << ", beta = " << beta )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dspmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dspmv( F77_UL, &F77_N, &alpha, AP, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dspmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dspmv( F77_UL, &F77_N, &alpha, AP, X, &F77_incX, &beta, Y, &F77_incY );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dspmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** spr */

template<>
void OpenMPBLAS2::spr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const float alpha,
    const float* X,
    const IndexType incX,
    float* Ap,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "spr<float>: N = " << N << ", alpha = " << alpha )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sspr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_sspr( F77_UL, &F77_N, &alpha, X, &F77_incX, Ap );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sspr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_sspr( F77_UL, &F77_N, &alpha, X, &F77_incX, Ap );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sspr", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::spr(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const double alpha,
    const double* X,
    const IndexType incX,
    double* Ap,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "spr<float>: N = " << N << ", alpha = " << alpha )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dspr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dspr( F77_UL, &F77_N, &alpha, X, &F77_incX, Ap );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dspr", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dspr( F77_UL, &F77_N, &alpha, X, &F77_incX, Ap );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dspr", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** spr2 */

template<>
void OpenMPBLAS2::spr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const float alpha,
    const float* X,
    const IndexType incX,
    const float* Y,
    const IndexType incY,
    float* A,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "spr2<float>: N = " << N << ", alpha = " << alpha )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sspr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_sspr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sspr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_sspr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sspr2", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
}

template<>
void OpenMPBLAS2::spr2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const IndexType N,
    const double alpha,
    const double* X,
    const IndexType incX,
    const double* Y,
    const IndexType incY,
    double* A,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "spr2<double>: N = " << N << ", alpha = " << alpha )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL;
#ifdef F77_CHAR
    F77_CHAR F77_UL;
#else
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX, F77_incY = incY;
#else
#define F77_N N
#define F77_incX incX
#define F77_incY incY
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dspr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dspr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dspr2", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
#endif
        F77_dspr2( F77_UL, &F77_N, &alpha, X, &F77_incX, Y, &F77_incY, A );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dspr2", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** tpmv */

template<>
void OpenMPBLAS2::tpmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const float* Ap,
    float* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tpmv<float>: N = " << N )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stpmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stpmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stpmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stpmv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stpmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stpmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stpmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stpmv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_stpmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::tpmv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const double* Ap,
    double* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tpmv<double>: N = " << N )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtpmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtpmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtpmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtpmv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtpmv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtpmv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtpmv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtpmv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtpmv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** tpsv */

template<>
void OpenMPBLAS2::tpsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const float* Ap,
    float* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tpsv<float>: N = " << N )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stpsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stpsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stpsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stpsv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_stpsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_stpsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_stpsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_stpsv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_stpsv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS2::tpsv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType N,
    const double* Ap,
    double* X,
    const IndexType incX,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger, "tpsv<double>: N = " << N )

    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA;
    char UL;
    char DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_incX = incX;
#else
#define F77_N N
#define F77_incX incX
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Uplo == CblasUpper )
        {
            UL = 'U';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtpsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtpsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtpsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtpsv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Uplo == CblasUpper )
        {
            UL = 'L';
        }
        else if ( Uplo == CblasLower )
        {
            UL = 'U';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtpsv", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasNoTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasTrans )
        {
            TA = 'N';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtpsv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( Diag == CblasUnit )
        {
            DI = 'U';
        }
        else if ( Diag == CblasNonUnit )
        {
            DI = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtpsv", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtpsv( F77_UL, F77_TA, F77_DI, &F77_N, Ap, X, &F77_incX );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtpsv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS2::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, gemv, double )

    // all other routines are not used in LAMA yet
}


/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS2 routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS2::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS2::initialized = registerInterface();

} /* namespace lama */

