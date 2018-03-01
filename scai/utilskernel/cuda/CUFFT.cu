/**
 * @file blaskernel/cuda/CUDABLAS1.cu
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @brief Wrapper implementations for BLAS1 routines in CUDA using cuBLAS
 * @author Lauretta Schubert, Thomas Brandes, Eric Stricker
 * @date 05.07.2012
 */

// hpp
#include <scai/utilskernel/cuda/CUFFT.hpp>

// local library
#include <scai/utilskernel/cuda/CUFFTWrapper.hpp>
#include <scai/utilskernel/FFTKernelTrait.hpp>
#include <scai/utilskernel/cuda/CUDAUtils.hpp>
#include <scai/utilskernel/cuda/CUDASparseUtils.hpp>

// internal scai libraries
#include <scai/tasking/cuda/CUDAStreamSyncToken.hpp>
#include <scai/kregistry/KernelRegistry.hpp>
#include <scai/hmemo.hpp>

#include <scai/common/BinaryOp.hpp>
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/cuda/CUDAAccess.hpp>
#include <scai/common/cuda/launchHelper.hpp>
#include <scai/common/macros/unused.hpp>
#include <scai/common/TypeTraits.hpp>

#include <cufft.h>

using namespace scai::tasking;
using scai::common::TypeTraits;

namespace scai
{

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( CUFFT::logger, "CUDA.utilskernel" )

/* ---------------------------- paddedFoward1D ---------------------------- */

template<typename ValueType>
void CUFFT::paddedForward1D(
    const IndexType n,
    const IndexType npad,
    const ValueType in[],
    common::Complex<ValueType> out[] )
{
    typedef common::Complex<ValueType> ComplexType;

    SCAI_CHECK_CUDA_ACCESS

    size_t size = sizeof(ComplexType) * npad;

    ComplexType* d_in;

    SCAI_CUDA_DRV_CALL( cuMemAlloc( ( CUdeviceptr* ) &d_in, size ), "cuMemAlloc failed" );

    ComplexType zero(0);

    CUDAUtils::setVal( d_in, npad, zero, common::BinaryOp::COPY );
    CUDASparseUtils::set( d_in, in, n, common::BinaryOp::COPY );

    cufftHandle plan;
    SCAI_CUFFT_CALL( cufftPlan1d(&plan, npad, CUFFTWrapper<ValueType>::getTypeC2C(), 1 ), "creation of plan for 1d fft failed")

    CUFFTWrapper<ValueType>::execute( plan, d_in, out, CUFFT_FORWARD );

    SCAI_CUDA_RT_CALL( cudaDeviceSynchronize(), "device synchronize")

    SCAI_CUFFT_CALL( cufftDestroy(plan), "destruction of plan failed" )

    SCAI_CUDA_DRV_CALL( cuMemFree( ( CUdeviceptr ) d_in ), "cuMemFree failed" )
}

/* ---------------------------- paddedBackwardD ---------------------------- */

template<typename ValueType>
void CUFFT::paddedBackward1D(
    const IndexType n,
    const IndexType npad,
    const ValueType in[],
    common::Complex<ValueType> out[] )
{
    typedef common::Complex<ValueType> ComplexType;

    SCAI_CHECK_CUDA_ACCESS

    size_t size = sizeof( ComplexType ) * npad;

    ComplexType* d_in; // static_cast<ComplexType*>( mem->allocate( size ) );

    SCAI_CUDA_DRV_CALL( cuMemAlloc( ( CUdeviceptr* ) &d_in, size ), "cuMemAlloc failed" );

    ComplexType zero(0);

    CUDAUtils::setVal( d_in, npad, zero, common::BinaryOp::COPY );
    CUDASparseUtils::set( d_in, in, n, common::BinaryOp::COPY );

    cufftHandle plan;
    SCAI_CUFFT_CALL( cufftPlan1d(&plan, npad, CUFFTWrapper<ValueType>::getTypeC2C(), 1 ), "creation of plan for 1d fft failed")

    CUFFTWrapper<ValueType>::execute( plan, d_in, out, CUFFT_INVERSE );

    SCAI_CUDA_RT_CALL( cudaDeviceSynchronize(), "device synchronize")

    SCAI_CUFFT_CALL( cufftDestroy(plan), "destruction of plan failed" )

    SCAI_CUDA_DRV_CALL( cuMemFree( ( CUdeviceptr ) d_in ), "cuMemFree failed" )
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void CUFFT::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    using kregistry::KernelRegistry;
    const common::ContextType ctx = common::ContextType::CUDA;
    SCAI_LOG_INFO( logger, "register FFT routines implemented by CUFFT in KernelRegistry [" << flag << "]" )
    KernelRegistry::set<FFTKernelTrait::paddedForward1D<ValueType> >( CUFFT::paddedForward1D, ctx, flag );
    KernelRegistry::set<FFTKernelTrait::paddedBackward1D<ValueType> >( CUFFT::paddedBackward1D, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUFFT::CUFFT()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_TYPELIST( float, double )>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUFFT::~CUFFT()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_TYPELIST( float, double )>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUFFT CUFFT::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
