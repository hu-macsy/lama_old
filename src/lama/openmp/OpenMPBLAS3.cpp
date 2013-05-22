/**
 * @file OpenMPBLAS3.cpp
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
 * @brief OpenMPBLAS3.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * $Id$
 */

// hpp
#include <lama/openmp/OpenMPBLAS3.hpp>

// others
#include <lama/BLASInterface.hpp>
#include <lama/LAMAInterfaceRegistry.hpp>

// macros
#include <lama/macros/unused.hpp>

namespace lama
{

LAMA_LOG_DEF_LOGGER( OpenMPBLAS3::logger, "OpenMP.BLAS3" )

template<>
void OpenMPBLAS3::gemm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_TRANSPOSE TransB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* B,
    const IndexType ldb,
    const float beta,
    float* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA, TB;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_TB;
#else
#define F77_TA &TA
#define F77_TB &TB
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = n, F77_K = k, F77_lda = lda, F77_ldb = ldb;
    F77_INT F77_ldc = ldc;
#else
#define F77_M m
#define F77_N n
#define F77_K k
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( TransB == CblasTrans )
        {
            TB = 'T';
        }
        else if ( TransB == CblasConjTrans )
        {
            TB = 'C';
        }
        else if ( TransB == CblasNoTrans )
        {
            TB = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_sgemm", "Illegal TransB setting, %d\n", TransB );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
        F77_TB = C2F_CHAR( &TB );
#endif
        F77_sgemm( F77_TA, F77_TB, &F77_M, &F77_N, &F77_K, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( TransA == CblasTrans )
        {
            TB = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TB = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TB = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( TransB == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransB == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransB == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
        F77_TB = C2F_CHAR( &TB );
#endif
        F77_sgemm( F77_TA, F77_TB, &F77_N, &F77_M, &F77_K, &alpha, B, &F77_ldb, A, &F77_lda, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sgemm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS3::gemm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_TRANSPOSE TransB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* B,
    const IndexType ldb,
    const double beta,
    double* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char TA, TB;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_TB;
#else
#define F77_TA &TA
#define F77_TB &TB
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = n, F77_K = k, F77_lda = lda, F77_ldb = ldb;
    F77_INT F77_ldc = ldc;
#else
#define F77_M m
#define F77_N n
#define F77_K k
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgemm", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( TransB == CblasTrans )
        {
            TB = 'T';
        }
        else if ( TransB == CblasConjTrans )
        {
            TB = 'C';
        }
        else if ( TransB == CblasNoTrans )
        {
            TB = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dgemm", "Illegal TransB setting, %d\n", TransB );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
        F77_TB = C2F_CHAR( &TB );
#endif
        F77_dgemm( F77_TA, F77_TB, &F77_M, &F77_N, &F77_K, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( TransA == CblasTrans )
        {
            TB = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TB = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TB = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgemm", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

        if ( TransB == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransB == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransB == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dgemm", "Illegal TransB setting, %d\n", TransB );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_TA = C2F_CHAR( &TA );
        F77_TB = C2F_CHAR( &TB );
#endif
        F77_dgemm( F77_TA, F77_TB, &F77_N, &F77_M, &F77_K, &alpha, B, &F77_ldb, A, &F77_lda, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dgemm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** symm */

template<>
void OpenMPBLAS3::symm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_SIDE Side,
    const enum CBLAS_UPLO Uplo,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* B,
    const IndexType ldb,
    const float beta,
    float* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char SD, UL;
#ifdef F77_CHAR
    F77_CHAR F77_SD, F77_UL;
#else
#define F77_SD &SD
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = n, F77_lda = lda, F77_ldb = ldb;
    F77_INT F77_ldc = ldc;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Side == CblasRight )
        {
            SD = 'R';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssymm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssymm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_SD = C2F_CHAR( &SD );
#endif
        F77_ssymm( F77_SD, F77_UL, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Side == CblasRight )
        {
            SD = 'L';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssymm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssymm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_SD = C2F_CHAR( &SD );
#endif
        F77_ssymm( F77_SD, F77_UL, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssymm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS3::symm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_SIDE Side,
    const enum CBLAS_UPLO Uplo,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* B,
    const IndexType ldb,
    const double beta,
    double* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char SD, UL;
#ifdef F77_CHAR
    F77_CHAR F77_SD, F77_UL;
#else
#define F77_SD &SD
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = N, F77_lda = lda, F77_ldb = ldb;
    F77_INT F77_ldc = ldc;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Side == CblasRight )
        {
            SD = 'R';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsymm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsymm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_SD = C2F_CHAR( &SD );
#endif
        F77_dsymm( F77_SD, F77_UL, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Side == CblasRight )
        {
            SD = 'L';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsymm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsymm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_SD = C2F_CHAR( &SD );
#endif
        F77_dsymm( F77_SD, F77_UL, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsymm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** trmm */

template<>
void OpenMPBLAS3::trmm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_SIDE Side,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* A,
    const IndexType lda,
    float* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TA, SD, DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_SD, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_SD &SD
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = n, F77_lda = lda, F77_ldb = ldb;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_ldb ldb
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Side == CblasRight )
        {
            SD = 'R';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strmm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strmm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strmm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_strmm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strmm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Side == CblasRight )
        {
            SD = 'L';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strmm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strmm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strmm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_strmm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strmm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_strmm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS3::trmm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_SIDE Side,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* A,
    const IndexType lda,
    double* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TA, SD, DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_SD, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_SD &SD
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = n, F77_lda = lda, F77_ldb = ldb;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_ldb ldb
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Side == CblasRight )
        {
            SD = 'R';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrmm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrmm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrmm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_dtrmm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrmm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Side == CblasRight )
        {
            SD = 'L';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrmm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrmm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrmm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_dtrmm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrmm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtrmm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** trsm */

template<>
void OpenMPBLAS3::trsm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_SIDE Side,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType m,
    const IndexType n,
    const float alpha,
    const float* A,
    const IndexType lda,
    float* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TA, SD, DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_SD, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_SD &SD
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_N = m, F77_N = n, F77_lda = lda, F77_ldb = ldb;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_ldb ldb
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Side == CblasRight )
        {
            SD = 'R';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strsm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_strsm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strsm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Side == CblasRight )
        {
            SD = 'L';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_strsm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_strsm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_strsm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_strsm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_strsm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_strsm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS3::trsm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_SIDE Side,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_DIAG Diag,
    const IndexType m,
    const IndexType n,
    const double alpha,
    const double* A,
    const IndexType lda,
    double* B,
    const IndexType ldb,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TA, SD, DI;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL, F77_SD, F77_DI;
#else
#define F77_TA &TA
#define F77_UL &UL
#define F77_SD &SD
#define F77_DI &DI
#endif
#ifdef F77_INT
    F77_INT F77_M = m, F77_N = n, F77_lda = lda, F77_ldb = ldb;
#else
#define F77_M m
#define F77_N n
#define F77_lda lda
#define F77_ldb ldb
#endif
    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if ( order == CblasColMajor )
    {
        if ( Side == CblasRight )
        {
            SD = 'R';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'L';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrsm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrsm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrsm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_dtrsm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        F77_dtrsm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_M, &F77_N, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else if ( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if ( Side == CblasRight )
        {
            SD = 'L';
        }
        else if ( Side == CblasLeft )
        {
            SD = 'R';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dtrsm", "Illegal Side setting, %d\n", Side );
            RowMajorStrg = 0;
            return;
        }

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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dtrsm", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( TransA == CblasTrans )
        {
            TA = 'T';
        }
        else if ( TransA == CblasConjTrans )
        {
            TA = 'C';
        }
        else if ( TransA == CblasNoTrans )
        {
            TA = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 4, "cblas_dtrsm", "Illegal Trans setting, %d\n", TransA );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 5, "cblas_dtrsm", "Illegal Diag setting, %d\n", Diag );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TA = C2F_CHAR( &TA );
        F77_SD = C2F_CHAR( &SD );
        F77_DI = C2F_CHAR( &DI );
#endif
        //          dtrsm(F77_SD, F77_UL, F77_TA, F77_DI, &F77_N, &F77_M, &alpha, A,
        //                                &F77_lda, B, &F77_ldb);
        F77_dtrsm( F77_SD, F77_UL, F77_TA, F77_DI, &F77_N, &F77_M, &alpha, A, &F77_lda, B, &F77_ldb );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dtrsm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** syrk */

template<>
void OpenMPBLAS3::syrk(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE Trans,
    const IndexType n,
    const IndexType k,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float beta,
    float* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TR;
#ifdef F77_CHAR
    F77_CHAR F77_TR, F77_UL;
#else
#define F77_TR &TR
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = N, F77_K = K, F77_lda = lda;
    F77_INT F77_ldc = ldc;
#else
#define F77_N n
#define F77_K k
#define F77_lda lda
#define F77_ldc ldc
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssyrk", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'T';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'C';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssyrk", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_ssyrk( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, &beta, C, &F77_ldc );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssyrk", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'T';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssyrk", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_ssyrk( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssyrk", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS3::syrk(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE Trans,
    const IndexType n,
    const IndexType k,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double beta,
    double* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TR;
#ifdef F77_CHAR
    F77_CHAR F77_TR, F77_UL;
#else
#define F77_TR &TR
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = n, F77_K = k, F77_lda = lda;
    F77_INT F77_ldc = ldc;
#else
#define F77_N n
#define F77_K k
#define F77_lda lda
#define F77_ldc ldc
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsyrk", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'T';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'C';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsyrk", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_dsyrk( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, &beta, C, &F77_ldc );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsyrk", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'T';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsyrk", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_dsyrk( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsyrk", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/** syrk2 */

template<>
void OpenMPBLAS3::syrk2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE Trans,
    const IndexType n,
    const IndexType k,
    const float alpha,
    const float* A,
    const IndexType lda,
    const float* B,
    const IndexType ldb,
    const float beta,
    float* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TR;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL;
#else
#define F77_TR &TR
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = n, F77_K = k, F77_lda = lda, F77_ldb = ldb;
    F77_INT F77_ldc = ldc;
#else
#define F77_N n
#define F77_K k
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_ssyr2k", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'T';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'C';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssyr2k", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_ssyr2k( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssyr2k", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'T';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_ssyr2k", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_ssyr2k( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_ssyr2k", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

template<>
void OpenMPBLAS3::syrk2(
    const enum CBLAS_ORDER order,
    const enum CBLAS_UPLO Uplo,
    const enum CBLAS_TRANSPOSE Trans,
    const IndexType n,
    const IndexType k,
    const double alpha,
    const double* A,
    const IndexType lda,
    const double* B,
    const IndexType ldb,
    const double beta,
    double* C,
    const IndexType ldc,
    SyncToken* syncToken )
{
    if ( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    char UL, TR;
#ifdef F77_CHAR
    F77_CHAR F77_TA, F77_UL;
#else
#define F77_TR &TR
#define F77_UL &UL
#endif
#ifdef F77_INT
    F77_INT F77_N = n, F77_K = k, F77_lda = lda, F77_ldb = ldb;
    F77_INT F77_ldc = ldc;
#else
#define F77_N n
#define F77_K k
#define F77_lda lda
#define F77_ldb ldb
#define F77_ldc ldc
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_dsyr2k", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'T';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'C';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'N';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsyr2k", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_dsyr2k( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
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
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsyr2k", "Illegal Uplo setting, %d\n", Uplo );
            RowMajorStrg = 0;
            return;
        }

        if ( Trans == CblasTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasConjTrans )
        {
            TR = 'N';
        }
        else if ( Trans == CblasNoTrans )
        {
            TR = 'T';
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 3, "cblas_dsyr2k", "Illegal Trans setting, %d\n", Trans );
            RowMajorStrg = 0;
            return;
        }

#ifdef F77_CHAR
        F77_UL = C2F_CHAR( &UL );
        F77_TR = C2F_CHAR( &TR );
#endif
        F77_dsyr2k( F77_UL, F77_TR, &F77_N, &F77_K, &alpha, A, &F77_lda, B, &F77_ldb, &beta, C, &F77_ldc );
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_dsyr2k", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS3::setInterface( BLASInterface& BLAS )
{
    // Note: macro takes advantage of same name for routines and type definitions 
    //       ( e.g. routine CUDABLAS1::sum<T> is set for BLAS::BLAS1::sum variable

    LAMA_INTERFACE_REGISTER_T( BLAS, gemm, float )
    LAMA_INTERFACE_REGISTER_T( BLAS, gemm, double )

    // trsm routines are not used yet by LAMA
}

/* --------------------------------------------------------------------------- */
/*    Static registration of the BLAS3 routines                                */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS3::registerInterface()
{
    LAMAInterface& interface = LAMAInterfaceRegistry::getRegistry().modifyInterface( Context::Host );
    setInterface( interface.BLAS );
    return true;
}

/* --------------------------------------------------------------------------- */
/*    Static initialiazion at program start                                    */
/* --------------------------------------------------------------------------- */

bool OpenMPBLAS3::initialized = registerInterface();

} /* namespace lama */
