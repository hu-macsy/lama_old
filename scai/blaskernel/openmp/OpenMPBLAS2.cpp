/**
 * @file OpenMPBLAS2.cpp
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
 * @brief Own implementation of BLAS2 kernels for Host using OpenMP parallelization.
 * @author Eric Schricker
 * @date 09.10.2013
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPBLAS2.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/tracing.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/bind.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

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
    SCAI_REGION( "OpenMP.BLAS2.gemv" )

    SCAI_LOG_INFO( logger,
                   "gemv<" << common::TypeTraits<ValueType>::id() << ">: M = " << M << ", N = " << N
                   << ", LDA = " << lda << ", incX = " << incX << ", incY = " << incY
                   << ", alpha = " << alpha << ", beta = " << beta )

    if ( M == 0 )
    {
        return; // empty X, Y, A
    }

    // N == 0: empty A, but deal with X, Y, we can handle this here
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "no asynchronous execution for openmp possible at this level." )
    }

    if ( order == CblasColMajor )
    {
        if ( TransA == CblasNoTrans )
        {
            //'N'
            // y = alpha * A * x + beta * y
            ValueType Z;

            if ( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) 

                for ( IndexType i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < N; j++ )
                    {
                        Z += A[lda * j + i] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) 
                for ( IndexType i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < N; j++ )
                    {
                        Z += A[lda * j + i] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }
        }
        else if ( TransA == CblasTrans )
        {
            //'T'
            // y = alpha * A^T * x + beta * y
            ValueType Z;

            if ( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) 

                for ( IndexType i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < M; j++ )
                    {
                        Z += A[lda * i + j] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) 
                for ( IndexType i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < M; j++ )
                    {
                        Z += A[lda * i + j] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }
        }
        else if ( TransA == CblasConjTrans )
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < N; i++ )
            {
                ValueType Z = 0;

                for ( IndexType j = 0; j < M; j++ )
                {
                    Z += common::Math::conj( A[lda * i + j] ) * X[j * incX];
                }

                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
        }
        else
        {
            COMMON_THROWEXCEPTION( "Illegal TransA setting " << TransA )
        }
    }
    else if ( order == CblasRowMajor )
    {
        if ( TransA == CblasNoTrans )
        {
            //'T'
            ValueType Z;

            if ( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) 

                for ( IndexType i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < N; j++ )
                    {
                        Z += A[lda * i + j] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) 
                for ( IndexType i = 0; i < M; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < N; j++ )
                    {
                        Z += A[lda * i + j] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }
        }
        else if ( TransA == CblasTrans )
        {
            //'N'
            ValueType Z;

            if ( incX == 1 && incY == 1 )
            {
                #pragma omp parallel for private(Z) 

                for ( IndexType i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < M; j++ )
                    {
                        Z += A[lda * j + i] * X[j];
                    }

                    Y[i] = Z * alpha + Y[i] * beta;
                }
            }
            else
            {
                //incX != 1 || incY != 1
                #pragma omp parallel for private(Z) 
                for ( IndexType i = 0; i < N; i++ )
                {
                    Z = static_cast<ValueType>( 0.0 );

                    for ( IndexType j = 0; j < M; j++ )
                    {
                        Z += A[lda * j + i] * X[j * incX];
                    }

                    Y[i * incY] = Z * alpha + Y[i * incY] * beta;
                }
            }
        }
        else if ( TransA == CblasConjTrans )
        {
            #pragma omp parallel for 

            for ( IndexType i = 0; i < N; i++ )
            {
                ValueType Z = 0;

                for ( IndexType j = 0; j < M; j++ )
                {
                    Z += common::Math::conj( A[lda * j + i] ) * X[j * incX];
                }

                Y[i * incY] = Z * alpha + Y[i * incY] * beta;
            }
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

template<typename ValueType>
void OpenMPBLAS2::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::Host;
    SCAI_LOG_DEBUG( logger, "register BLAS2 routines for OpenMP in Kernel Registry" )
    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( OpenMPBLAS2::gemv, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS2::OpenMPBLAS2()
{
    SCAI_LOG_INFO( logger, "register BLAS2 routines for OpenMP in Kernel Registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPBLAS2::~OpenMPBLAS2()
{
    SCAI_LOG_INFO( logger, "unregister BLAS2 routines for OpenMP in Kernel Registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS2 OpenMPBLAS2::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
