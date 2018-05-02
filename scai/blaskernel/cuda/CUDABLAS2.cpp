/**
 * @file CUDABLAS2.cpp
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
 * @brief CUDA implementations of BLAS2 routines for the class CUDABLAS2.
 * @author Lauretta Schubert
 * @date 05.07.2012
 */

// hpp
#include <scai/blaskernel/cuda/CUDABLAS2.hpp>

// local library
#include <scai/blaskernel/cuda/CUBLASWrapper.hpp>
#include <scai/blaskernel/cuda/CUBLASTrait.hpp>
#include <scai/blaskernel/BLASKernelTrait.hpp>

// internal scai libraries
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>

using namespace scai::tasking;

namespace scai
{

using common::TypeTraits;

namespace blaskernel
{

SCAI_LOG_DEF_LOGGER( CUDABLAS2::logger, "CUDA.BLAS2" )

/* ---------------------------------------------------------------------------------------*/
/*    gemv                                                                                */
/* ---------------------------------------------------------------------------------------*/

/** gemv */

template<typename ValueType>
void CUDABLAS2::gemv(
    const common::MatrixOp op,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* const A,
    const IndexType lda,
    const ValueType* const x,
    const IndexType incx,
    const ValueType beta,
    ValueType* const y,
    const IndexType incy )
{
    // Attention: cuBLAS expects the matrix A stored in column-major format, 
    //            so we deal the transposed problem here 

    CUBLASTrait::BLASTrans transA = CUBLASTrait::castTrans( common::combine( op, common::MatrixOp::TRANSPOSE ) );

    SCAI_CHECK_CUDA_ACCESS

    cudaStream_t stream = 0; // default stream if no syncToken is given
    CUDAStreamSyncToken* syncToken = CUDAStreamSyncToken::getCurrentSyncToken();

    if ( syncToken )
    {
        stream = syncToken->getCUDAStream();
    }

    cublasHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuBLASHandle();

    SCAI_CUBLAS_CALL( cublasSetStream( handle, stream ),
                      "CUDABLAS2::gemv set cublas kernel stream = " << stream );
    SCAI_LOG_INFO( logger,
                   "gemv<" << TypeTraits<ValueType>::id() << "> with cuBLAS: m = " << m << " x " << n )

    // as we work on transposed data we swap here m and n

    CUBLASWrapper<ValueType>::gemv( handle, transA,  n,  m, alpha, A,  lda, x, incx , beta, y, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }

    SCAI_CUBLAS_CALL( cublasSetStream( handle, NULL ), "CUDABLAS2::gemv set stream NULL" );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDABLAS2::geam(
    ValueType* C,
    const IndexType ldc,
    const IndexType m,
    const IndexType n,
    const ValueType alpha,
    const ValueType* A,
    const IndexType lda,
    const common::MatrixOp opA,
    const ValueType beta,
    const ValueType* B,
    const IndexType ldb,
    const common::MatrixOp opB )
{
    cublasHandle_t handle = common::CUDAAccess::getCurrentCUDACtx().getcuBLASHandle();

    typedef CUBLASTrait::BLASTrans BLASTrans;

    BLASTrans transA = CUBLASTrait::castTrans( opA );
    BLASTrans transB = CUBLASTrait::castTrans( opB );

    SCAI_CHECK_CUDA_ACCESS

    if ( m == 0 || n == 0 )
    {
        return;    // important as call of geam might return an error 
    }

    cudaStream_t stream = 0; // default stream if no syncToken is given

    SCAI_CUBLAS_CALL( cublasSetStream( handle, stream ),
                      "CUDABLAS2::geam set cublas kernel stream = " << stream );

    SCAI_LOG_INFO( logger,
                    "geam<" << TypeTraits<ValueType>::id() << "> with cuBLAS: "
                     << "C ( " << m << " x " << n << ", ldc = " << ldc << " ) "
                     << " = " << alpha << " * A ( lda = " << lda << " ) " << opA 
                     << " + " << beta << " * B ( ldb = " << ldb << " ) " << opB )

    // we swap m and n as geam expects matrices in columns-major format

    CUBLASWrapper<ValueType>::geam( handle, transA, transB, n, m,
                                    alpha, A, lda, 
                                    beta, B, ldb,
                                    C, ldc );

    SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDABLAS2::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register BLAS2 routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )
    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( CUDABLAS2::gemv, ctx, flag );
    KernelRegistry::set<BLASKernelTrait::geam<ValueType> >( CUDABLAS2::geam, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDABLAS2::CUDABLAS2()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUDABLAS2::~CUDABLAS2()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_CUDA_LIST>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUDABLAS2 CUDABLAS2::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
