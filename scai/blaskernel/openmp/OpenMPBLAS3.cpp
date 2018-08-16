/**
 * @file OpenMPBLAS3.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Own implementation of BLAS3 kernels for Host using OpenMP parallelization.
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
#include <scai/tracing.hpp>

namespace scai
{

using tasking::TaskSyncToken;
using common::TypeTraits;
using common::MatrixOp;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPBLAS3::logger, "OpenMP.BLAS3" )

template<typename ValueType>
void OpenMPBLAS3::gemm(
    const MatrixOp opA,
    const MatrixOp opB,
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
    SCAI_REGION( "OpenMP.BLAS3.gemm" )

    TaskSyncToken* syncToken = TaskSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        SCAI_LOG_WARN( logger, "asynchronous execution not supported yet" )
    }

    SCAI_LOG_INFO( logger,
                    "gemm<" << TypeTraits<ValueType>::id() << ">: "
                    << " C = " << alpha << " * A * B + " << beta << " * C"
                    << ", C is " << m << " x " << n << ", ldc = " << ldc 
                    << ", A is " << m << " x " << k << ", lda = " << lda << ", opA = " << opA
                    << ", B is " << k << " x " << n << ", ldb = " << ldb << ", opB = " << opB )

    if ( opA == MatrixOp::TRANSPOSE )
    {
        if ( opB == MatrixOp::NORMAL )
        {
            #pragma omp parallel for collapse(2)

            for ( IndexType h = 0; h < n; h++ )
            {
                for ( IndexType i = 0; i < m; i++ )
                {
                    ValueType temp = 0;

                    for ( IndexType j = 0; j < k; j++ )
                    {
                        temp += A[lda * j + i] * B[ldb * j + h];
                    }

                    C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                }
            }
        }
        else if ( opB == MatrixOp::CONJ_TRANSPOSE )
        {
            COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
        }
        else if ( opB == MatrixOp::TRANSPOSE )
        {
            #pragma omp parallel for collapse(2)

            for ( IndexType h = 0; h < n; h++ )
            {
                for ( IndexType i = 0; i < m; i++ )
                {
                    ValueType temp = 0;

                    for ( IndexType j = 0; j < k; j++ )
                    {
                        temp += A[lda * j + i] * B[ldb * h + j];
                    }

                    C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                }
            }
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal opA setting " << opA )
        }
    }
    else if ( opA == MatrixOp::NORMAL )
    {
        if ( opB == MatrixOp::NORMAL )
        {
            // A = 'N'; B = 'N'
            //std::cout << "lda:" << lda << ", ldb:" << ldb << ", ldc:" << ldc << "\n";
            //std::cout << "n:" << n << ", m:" << m << ", k:" << k << "\n";
            #pragma omp parallel for collapse(2) 

            for ( IndexType h = 0; h < n; h++ )
            {
                for ( IndexType i = 0; i < m; i++ )
                {
                    ValueType temp = 0;

                    for ( IndexType j = 0; j < k; j++ )
                    {
                        temp += A[lda * i + j] * B[ldb * j + h];
                    } 

                    C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                }
            }
        }
        else if ( opB == MatrixOp::TRANSPOSE )
        {
            #pragma omp parallel for collapse(2)

            for ( IndexType h = 0; h < n; h++ )
            {
                for ( IndexType i = 0; i < m; i++ )
                {
                    ValueType temp = 0;

                    for ( IndexType j = 0; j < k; j++ )
                    {
                        temp += A[lda * i + j] * B[ldb * h + j];
                    }

                    C[ldc * i + h] = alpha * temp + beta * C[ldc * i + h];
                }
            }
        }
        else if ( opB == MatrixOp::CONJ_TRANSPOSE )
        {
            COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal opA setting " << opA )
        }
    }
    else if ( opA == MatrixOp::CONJ_TRANSPOSE )
    {
        // Todo: implement
        COMMON_THROWEXCEPTION( "gemm for complexe matrix is not supported yet" )
        /*
        if( opB == MatrixOp::NORMAL )
        {

        }
        else if( opB == MatrixOp::CONJ_TRANSPOSE )
        {

        }
        else if( opB == MatrixOp::TRANSPOSE )
        {

        }
        else
        {
            COMMON_THROWEXCEPTION( "illegal transA setting " << TransA )
        }*/
    }

    return;
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPBLAS3::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::Host;
    SCAI_LOG_DEBUG( logger, "set BLAS3 routines for OpenMP in Interface" )
    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( OpenMPBLAS3::gemm, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPBLAS3::OpenMPBLAS3()
{
    SCAI_LOG_INFO( logger, "register BLAS3 routines for OpenMP in Kernel Registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

OpenMPBLAS3::~OpenMPBLAS3()
{
    SCAI_LOG_INFO( logger, "unregister BLAS3 routines for OpenMP in Kernel Registry" )

    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPBLAS3 OpenMPBLAS3::guard;

} /* end namespace blaskernel */

} /* end namespace scai */
