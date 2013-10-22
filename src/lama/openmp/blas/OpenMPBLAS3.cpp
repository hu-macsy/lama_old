/**
 * @file OpenMPBLAS3.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @author Eric Schricker
 * @date 15.10.2013
 * @since 1.0.0
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

template<typename T>
void OpenMPBLAS3::gemm(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const enum CBLAS_TRANSPOSE TransB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const T alpha,
    const T* A,
    const IndexType UNUSED(lda),
    const T* B,
    const IndexType UNUSED(ldb),
    const T beta,
    T* C,
    const IndexType UNUSED(ldc),
    SyncToken* syncToken )
{
    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if( order == CblasColMajor )
    {
        if( TransA == CblasTrans )
        {
            //'T'
            if( TransB == CblasNoTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[k * h + j] * B[k * i + j];
                        }
                        C[m * i + h] = alpha * temp + beta * C[m * i + h];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                LAMA_THROWEXCEPTION( "gemm for complexe matrix is not supported yet")
            }
            else if( TransB == CblasTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[k * h + j] * B[n * j + i];
                        }
                        C[m * i + h] = alpha * temp + beta * C[m * i + h];
                    }
                }
            }
            else
            {
                BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
                RowMajorStrg = 0;
                return;
            }
        }
        else if( TransA == CblasConjTrans )
        {
            LAMA_THROWEXCEPTION( "gemm for complexe matrix is not supported yet")
            //'C'
            if( TransB == CblasNoTrans )
            {

            }
            else if( TransB == CblasConjTrans )
            {

            }
            else if( TransB == CblasTrans )
            {

            }
            else
            {
                BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
                RowMajorStrg = 0;
                return;
            }
        }
        else if( TransA == CblasNoTrans )
        {
            if( TransB == CblasNoTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[m * j + h] * B[k * i + j];
                        }
                        C[m * i + h] = alpha * temp + beta * C[m * i + h];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                LAMA_THROWEXCEPTION( "gemm for complexe matrix is not supported yet")
            }
            else if( TransB == CblasTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[m * j + h] * B[n * j + i];
                        }
                        C[m * i + h] = alpha * temp + beta * C[m * i + h];
                    }
                }
            }
            else
            {
                BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
                RowMajorStrg = 0;
                return;
            }
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }
    }
    else if( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if( TransA == CblasTrans )
        {
            if( TransB == CblasNoTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[m * j + h] * B[n * j + i];
                        }
                        C[n * h + i] = alpha * temp + beta * C[n * h + i];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                LAMA_THROWEXCEPTION( "gemm for complexe matrix is not supported yet")
            }
            else if( TransB == CblasTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[m * j + h] * B[k * i + j];
                        }
                        C[n * h + i] = alpha * temp + beta * C[n * h + i];
                    }
                }
            }
            else
            {
                BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
                RowMajorStrg = 0;
                return;
            }
        }
        else if( TransA == CblasNoTrans )
        {
            if( TransB == CblasNoTrans )
            {
                //A = 'N'; B = 'N'
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[k * h + j] * B[n * j + i];
                        }
                        C[n * h + i] = alpha * temp + beta * C[n * h + i];
                    }
                }
            }
            else if( TransB == CblasTrans )
            {
                T temp = 0.0;
#pragma omp parallel for collapse(2) private(temp)
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = 0.0;
                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[k * h + j] * B[k * i + j];
                        }
                        C[n * h + i] = alpha * temp + beta * C[n * h + i];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                LAMA_THROWEXCEPTION( "gemm for complexe matrix is not supported yet")
            }
            else
            {
                BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
                RowMajorStrg = 0;
                return;
            }
        }
        else if( TransA == CblasConjTrans )
        {
            LAMA_THROWEXCEPTION( "gemm for complexe matrix is not supported yet")
            if( TransB == CblasNoTrans )
            {

            }
            else if( TransB == CblasConjTrans )
            {

            }
            else if( TransB == CblasTrans )
            {

            }
            else
            {
                BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemm", "Illegal TransA setting, %d\n", TransA );
                RowMajorStrg = 0;
                return;
            }
        }
    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sgemm", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS3::setInterface( BLASInterface& BLAS )
{
    LAMA_LOG_INFO( logger, "set BLAS3 routines for OpenMP in Interface" )

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
