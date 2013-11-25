/**
 * @file OpenMPBLAS2.cpp
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
 * @brief OpenMPBLAS2.cpp
 * @author Eric Schricker
 * @date 09.10.2013
 * @since 1.0.0
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

template<typename T>
void OpenMPBLAS2::gemv(
    const enum CBLAS_ORDER order,
    const enum CBLAS_TRANSPOSE TransA,
    const IndexType M,
    const IndexType N,
    const T alpha,
    const T* A,
    const IndexType lda,
    const T* X,
    const IndexType incX,
    const T beta,
    T* Y,
    const IndexType incY,
    SyncToken* syncToken )
{
    LAMA_LOG_INFO( logger,
                   "gemv<float>: M = " << M << ", N = " << N << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY << ", alpha = " << alpha << ", beta = " << beta )

    if( M == 0 )
    {
        return; // empty X, Y, A
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here

    if( syncToken )
    {
        LAMA_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }


    IndexType RowMajorStrg;
    RowMajorStrg = 0;

    if( order == CblasColMajor )
    {
        if( TransA == CblasNoTrans )
        {
            //'N'
            // y = alpha * A * x + beta * y
            T Z = 0.0;

#pragma omp parallel for private(Z)
            for( int i = 0; i < M; i++ )
            {
                Z = 0.0;
                for( int j = 0; j < N; j++ )
                {
                    Z += A[lda * j + i] * X[j * incX];
                }
                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
        else if( TransA == CblasTrans )
        {
            //'T'
            // y = alpha * A^T * x + beta * y
            T Z = 0.0;

#pragma omp parallel for private(Z)
            for( int i = 0; i < N; i++ )
            {
                Z = 0.0;
                for( int j = 0; j < M; j++ )
                {
                    Z += A[lda * i + j] * X[j * incX];
                }
                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
        else if( TransA == CblasConjTrans )
        {
            //'C'
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
        }

    }
    else if( order == CblasRowMajor )
    {
        RowMajorStrg = 1;

        if( TransA == CblasNoTrans )
        {
            //'T'
            T Z = 0.0;

#pragma omp parallel for private(Z)
            for( int i = 0; i < M; i++ )
            {
                Z = 0.0;
                for( int j = 0; j < N; j++ )
                {
                    Z += A[lda * i + j] * X[j * incX];
                }
                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
        else if( TransA == CblasTrans )
        {
            //'N'
            T Z = 0.0;

#pragma omp parallel for private(Z)
            for( int i = 0; i < N; i++ )
            {
                Z = 0.0;
                for( int j = 0; j < M; j++ )
                {
                    Z += A[lda * j + i] * X[j * incX];
                }
                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
        else if( TransA == CblasConjTrans )
        {
            //TA = 'N'
        }
        else
        {
            BLASHelper::XERBLA_cpu( RowMajorStrg, 2, "cblas_sgemv", "Illegal TransA setting, %d\n", TransA );
            RowMajorStrg = 0;
            return;
        }

    }
    else
    {
        BLASHelper::XERBLA_cpu( RowMajorStrg, 1, "cblas_sgemv", "Illegal order setting, %d\n", order );
    }

    RowMajorStrg = 0;
    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS2::setInterface( BLASInterface& BLAS )
{
    LAMA_LOG_INFO( logger, "set BLAS2 routines for OpenMP in Interface" )

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

