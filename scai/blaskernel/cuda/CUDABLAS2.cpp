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
    const CBLAS_ORDER order,
    const CBLAS_TRANSPOSE trans,
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
    typedef CUBLASTrait::BLASTrans BLASTrans;
    IndexType order_m = m;
    IndexType order_n = n;
//    char trans_char = ' ';
    BLASTrans trans_char;

    //switch stuff because columnmajor to rowmajor
    if ( order == CblasRowMajor )
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = CUBLAS_OP_T;
        }
        else
        {
            trans_char = CUBLAS_OP_N;
        }

        order_m = n;
        order_n = m;
    }
    else
    {
        if ( trans == CblasNoTrans )
        {
            trans_char = CUBLAS_OP_N;
        }
        else
        {
            trans_char = CUBLAS_OP_T;
        }
    }

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
                   "gemv<" << TypeTraits<ValueType>::id() << "> with cuBLAS: m = " << order_m << " x " << order_n )
    CUBLASWrapper<ValueType>::gemv( handle, trans_char,  order_m ,  order_n , alpha, A,  lda, x, incx , beta, y, incy );

    // No error check here possible as kernel is started asynchronously

    if ( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }

    SCAI_CUBLAS_CALL( cublasSetStream( handle, NULL ), "CUDABLAS2::gemv set stream NULL" );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUDABLAS2::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::context::ContextType ctx = common::context::CUDA;
    SCAI_LOG_INFO( logger, "register BLAS2 routines implemented by CuBLAS in KernelRegistry [" << flag << "]" )
    KernelRegistry::set<BLASKernelTrait::gemv<ValueType> >( CUDABLAS2::gemv, ctx, flag );
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
