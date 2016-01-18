/**
 * @file OpenMPBLAS2.cpp
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
 * @brief Implementation of BLAS2 routines used in LAMA with own C++/OpenMP implementations.
 * @author Eric Schricker
 * @date 09.10.2013
 * @since 1.1.0
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPBLAS2.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/TypeTraits.hpp>

// boost
#include <boost/preprocessor.hpp>

namespace scai
{

using tasking::TaskSyncToken;
using common::TypeTraits;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPBLAS2::logger, "OpenMP.BLAS2" )

/** gemv */

template<typename ValueType>
void OpenMPBLAS2::gemv(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE TransA,
    const IndexType M,
    const IndexType N,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const ValueType* X,
    const IndexType incX,
    const ValueType beta,
    ValueType* Y,
    const IndexType incY )
{
    SCAI_LOG_INFO( logger,
                   "gemv<" << common::TypeTraits<ValueType>::id()<< ">: M = " << M << ", N = " << N 
                      << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY 
                      << ", alpha = " << alpha << ", beta = " << beta )

    if( M == 0 )
    {
        return; // empty X, Y, A
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if( order == CblasColMajor )
    {
        if( TransA == CblasNoTrans )
        {
            //'N'
            // y = alpha * A * x + beta * y
            ValueType Z;

            if( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )

                for( int i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < N; j++ )
                    {
                        Z += A[lda * j + i] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )
                for( int i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < N; j++ )
                    {
                        Z += A[lda * j + i] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }

        }
        else if( TransA == CblasTrans )
        {
            //'T'
            // y = alpha * A^T * x + beta * y
            ValueType Z;

            if( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )

                for( int i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < M; j++ )
                    {
                        Z += A[lda * i + j] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )
                for( int i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < M; j++ )
                    {
                        Z += A[lda * i + j] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }

        }
        else if( TransA == CblasConjTrans )
        {
            //'C'
        }
        else
        {
        	COMMON_THROWEXCEPTION( "Illegal TransA setting " << TransA )
        }

    }
    else if( order == CblasRowMajor )
    {
        if( TransA == CblasNoTrans )
        {
            //'T'
            ValueType Z;

            if( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )

                for( int i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < N; j++ )
                    {
                        Z += A[lda * i + j] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )
                for( int i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < N; j++ )
                    {
                        Z += A[lda * i + j] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }

        }
        else if( TransA == CblasTrans )
        {
            //'N'
            ValueType Z;

            if( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )

                for( int i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < M; j++ )
                    {
                        Z += A[lda * j + i] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) schedule( SCAI_OMP_SCHEDULE )
                for( int i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>(0.0);

                    for( int j = 0; j < M; j++ )
                    {
                        Z += A[lda * j + i] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }
        }
        if( TransA == CblasConjTrans )
        {
            //TA = 'N'
        }
        else
        {
        	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
        }

    }
    else
    {
    	COMMON_THROWEXCEPTION( "illegal order setting " << order )
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void OpenMPBLAS2::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;
    using common::context::Host;

    SCAI_LOG_INFO( logger, "register BLAS2 routines for OpenMP in Kernel Registry" )

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD ;   // lower priority

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

#define LAMA_BLAS2_REGISTER(z, I, _)                                                          \
    KernelRegistry::set<BLASKernelTrait::gemv<ARITHMETIC_HOST_TYPE_##I> >( gemv, Host, flag ); \

    BOOST_PP_REPEAT( ARITHMETIC_HOST_TYPE_CNT, LAMA_BLAS2_REGISTER, _ )

#undef LAMA_BLAS2_REGISTER

    // all other routines are not used in LAMA yet
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS2::OpenMPBLAS2()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

OpenMPBLAS2::~OpenMPBLAS2()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS2 OpenMPBLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
