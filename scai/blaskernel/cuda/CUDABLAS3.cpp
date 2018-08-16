/**
 * @file CUDABLAS3.cpp
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
 * @brief CUDA implementations of BLAS3 routines.
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/cuda/CUDABLAS3.hpp>

// local library
#include <scai/blaskernel/BLASKernelTrait.hpp>
#include <scai/blaskernel/cuda/CUBLASWrapper.hpp>

// internal scai library
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>

#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

#include <scai/tracing.hpp>

namespace scai
{

using namespace tasking;
using common::TypeTraits;
using common::MatrixOp;

/* ---------------------------------------------------------------------------------------*/
namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( CUDABLAS3::logger, "CUDA.BLAS3" )

/*    gemm                                                                                */
/* ---------------------------------------------------------------------------------------*/

template<typename ValueType>
void CUDABLAS3::gemm(
    const MatrixOp opA,
    const MatrixOp opB,
    const IndexType m,
    const IndexType n,
    const IndexType k,
    const ValueType alpha,
    const ValueType* const A,
    const IndexType lda,
    const ValueType* const B,
    const IndexType ldb,
    const ValueType beta,
    ValueType* const C,
    const IndexType ldc )
{
    SCAI_REGION( "CUDA.BLAS3.gemm" )

    SCAI_LOG_INFO( logger, "gemm<" << TypeTraits<ValueType>::id() << ">( m = " << m << ", n = " << n << ", k = " << k )

    CUBLASTrait::BLASTrans transA = CUBLASTrait::castTrans( opA );
    CUBLASTrait::BLASTrans transB = CUBLASTrait::castTrans( opB );

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    cublasHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuBLASHandle();

    SCAI_CUBLAS_CALL( cublasSetStream( handle, stream ),
                      "CUDABLAS3::gemm set cublas kernel stream = " << stream );

    SCAI_LOG_INFO( logger, "cublas gemm: m = " << m << " x " << n )

    // Attention: gemm of cuBLAS can only deal with column-major format so we
    //            compute transpose( C ) = alpha * transpose( op( B ) ) * transpose( op( A ) ) + beta * transpose( C )

    CUBLASWrapper<ValueType>::gemm( handle, transB, transA,  
                                    n, m,  k , alpha, B,  ldb , A,  lda , beta, C, ldc );

    // No error check here possible as kernel is started asynchronously in any case

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "stream = " << stream );
    }

    SCAI_CUBLAS_CALL( cublasSetStream( handle, NULL ), "CUDABLAS3::gemm set stream to NULL" );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDABLAS3::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register BLAS3 routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )
    KernelRegistry::set<BLASKernelTrait::gemm<ValueType> >( CUDABLAS3::gemm, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDABLAS3::CUDABLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUDABLAS3::~CUDABLAS3()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUDABLAS3 CUDABLAS3::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
