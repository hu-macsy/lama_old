/**
 * @file CUDABLAS2.cpp
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
 * @brief CUDABLAS2.cpp
 * @author Lauretta Schubert
 * @date 05.07.2012
 * @since 1.0.0
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
#include <scai/common/preprocessor.hpp>

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
    if( order == CblasRowMajor )
    {
        if( trans == CblasNoTrans )
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
        if( trans == CblasNoTrans )
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

    cublasHandle_t handle = common::CUDAAccess::getCurrentCUDADevice().getcuBLASHandle();

    SCAI_CUBLAS_CALL( cublasSetStream( handle, stream ),
                      "CUDABLAS2::gemv set cublas kernel stream = " << stream );

    SCAI_LOG_INFO( logger,
                   "gemv<" << TypeTraits<ValueType>::id() << "> with cuBLAS: m = " << order_m << " x " << order_n )

    CUBLASWrapper<ValueType>::gemv( handle, trans_char,  order_m ,  order_n , alpha, A,  lda, x, incx , beta, y, incy );

    // No error check here possible as kernel is started asynchronously

    if( !syncToken )
    {
        SCAI_CUDA_RT_CALL( cudaStreamSynchronize( stream ), "cudaStreamSynchronize( stream = " << stream << " )" );
    }

    SCAI_CUBLAS_CALL( cublasSetStream( handle, NULL ), "CUDABLAS2::gemv set stream NULL" );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

void CUDABLAS2::registerKernels( bool deleteFlag )
{
    using kregistry::KernelRegistry;

    const common::context::ContextType ctx = common::context::CUDA;

    KernelRegistry::KernelRegistryFlag flag = KernelRegistry::KERNEL_ADD;

    if ( deleteFlag )
    {
        flag = KernelRegistry::KERNEL_ERASE;
    }

    SCAI_LOG_INFO( logger, "register BLAS2 routines implemented by CuBLAS in KernelRegistry" )

    // register for one CUDA type: ARITHMETIC_CUDA_TYPE_xxx

#define LAMA_BLAS2_REGISTER(z, I, _)                                                            \
    KernelRegistry::set<BLASKernelTrait::gemv<ARITHMETIC_CUDA_TYPE_##I> >( gemv, ctx, flag );   \

    BOOST_PP_REPEAT( ARITHMETIC_CUDA_TYPE_CNT, LAMA_BLAS2_REGISTER, _ )

#undef LAMA_BLAS2_REGISTER
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUDABLAS2::CUDABLAS2()
{
    bool deleteFlag = false;
    registerKernels( deleteFlag );
}

CUDABLAS2::~CUDABLAS2()
{
    bool deleteFlag = true;
    registerKernels( deleteFlag );
}

CUDABLAS2 CUDABLAS2::guard;    // guard variable for registration

} /* end namespace blaskernel */

} /* end namespace scai */
