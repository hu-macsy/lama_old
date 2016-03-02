/**
 * @file OpenMPBLAS3.cpp
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
 * @brief Implementation of BLAS3 methods used in LAMA with own C++ implementation.
 * @author Eric Schricker
 * @date 15.10.2013
 * @since 1.1.0
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPBLAS3.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/mepr/Container.hpp>

namespace scai
{

using tasking::TaskSyncToken;
using common::TypeTraits;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPBLAS3::logger, "OpenMP.BLAS3" )

template<typename ValueType>
void OpenMPBLAS3::gemm(
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE TransA,
    const CBLAS_TRANSPOSE TransB,
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
    const IndexType ldc )
{
    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet" )
    }

    SCAI_LOG_INFO( logger,
                   "gemm<" << TypeTraits<ValueType>::id() << ">: " << "m = " << m << ", n = " << n 
                     << ", k = " << k << ", lda = " << lda << ", ldb = " << ldb << ", ldc = " << ldc 
                     << ", alpha = " << alpha << ", beta = " << beta )

    if( order == CblasColMajor )
    {
        if( TransA == CblasTrans )
        {
            //'T'
            if( TransB == CblasNoTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * h + j] * B[ldb * i + j];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
            }
            else if( TransB == CblasTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * h + j] * B[ldb * j + i];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else
            {
            	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
            }
        }
        else if( TransA == CblasConjTrans )
        {
            COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )

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
            	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
            }
        }
        else if( TransA == CblasNoTrans )
        {
            if( TransB == CblasNoTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * j + h] * B[ldb * i + j];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
            }
            else if( TransB == CblasTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * j + h] * B[ldb * j + i];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else
            {
            	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
                return;
            }
        }
        else
        {
        	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
            return;
        }
    }
    else if( order == CblasRowMajor )
    {
        if( TransA == CblasTrans )
        {
            if( TransB == CblasNoTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * j + i] * B[ldb * j + h];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
            }
            else if( TransB == CblasTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * j + i] * B[ldb * h + j];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else
            {
            	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
            }
        }
        else if( TransA == CblasNoTrans )
        {
            if( TransB == CblasNoTrans )
            {
                // A = 'N'; B = 'N'
                //std::cout << "lda:" << lda << ", ldb:" << ldb << ", ldc:" << ldc << "\n";
                //std::cout << "n:" << n << ", m:" << m << ", k:" << k << "\n";
                ValueType temp;

                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )
                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * i + j] * B[ldb * j + h];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else if( TransB == CblasTrans )
            {
                ValueType temp;
                #pragma omp parallel for collapse(2) private(temp) schedule( SCAI_OMP_SCHEDULE )

                for( int h = 0; h < n; h++ )
                {
                    for( int i = 0; i < m; i++ )
                    {
                        temp = static_cast<ValueType>(0.0);

                        for( int j = 0; j < k; j++ )
                        {
                            temp += A[lda * i + j] * B[ldb * h + j];
                        }

                        C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                    }
                }
            }
            else if( TransB == CblasConjTrans )
            {
                COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
            }
            else
            {
            	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
            }
        }
        else if( TransA == CblasConjTrans )
        {
            COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )

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
            	COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
            }
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
void OpenMPBLAS3::RegistratorV<ValueType>::initAndReg( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using common::context::Host;
    using kregistry::KernelRegistry;

    SCAI_LOG_INFO( logger, "register BLAS3 OpenMP-routines for Host at kernel registry [" << flag << "]" )

    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( OpenMPBLAS3::gemm, Host, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS3::OpenMPBLAS3()
{
    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_HOST> ValueTypes;

    common::mepr::instantiate( kregistry::KernelRegistry::KERNEL_ADD, ValueTypes() );
}

OpenMPBLAS3::~OpenMPBLAS3()
{
    typedef common::mepr::ContainerV<RegistratorV, ARITHMETIC_HOST> ValueTypes;

    common::mepr::instantiate( kregistry::KernelRegistry::KERNEL_ERASE, ValueTypes() );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS3 OpenMPBLAS3::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
