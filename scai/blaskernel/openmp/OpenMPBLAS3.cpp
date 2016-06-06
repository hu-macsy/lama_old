/**
 * @file OpenMPBLAS3.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of BLAS3 methods used in LAMA with own C++ implementation.
 * @author Eric Schricker
 * @date 15.10.2013
 */

// hpp
#include <scai/blaskernel/openmp/OpenMPBLAS3.hpp>

#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/TaskSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

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
            // Todo: implement
            COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )

            //'C'
            /*
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
            }*/
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
            }
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
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
            // Todo: implement
            COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )

            /*
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
            }*/
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
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::Host;

    SCAI_LOG_INFO( logger, "set BLAS3 routines for OpenMP in Interface" )

    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( OpenMPBLAS3::gemm, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS3::OpenMPBLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_HOST_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPBLAS3::~OpenMPBLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_ARITHMETIC_HOST_LIST>::call(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS3 OpenMPBLAS3::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
