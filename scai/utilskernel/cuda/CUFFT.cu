/**
 * @file utilskernel/cuda/CUFFT.cu
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

#include <scai/tracing.hpp>

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

/* ----------------------------  fft  -------------------------------------- */

template<typename ValueType>
void CUFFT::fft( common::Complex<ValueType> x[], IndexType nb, IndexType n, IndexType m, int dir )
{
    SCAI_REGION( "CUDA.fft" )

    SCAI_LOG_INFO( logger, "fft @ CUDA, " << nb << " x " << n << " = 2 ** " << m )

    cufftHandle plan;

    SCAI_CHECK_CUDA_ACCESS

    SCAI_CUFFT_CALL( cufftPlan1d( &plan, n, CUFFTWrapper<ValueType>::getTypeC2C(), nb ), 
                     "creation of plan for 1d fft failed")

    if ( dir == 1 )
    {
        CUFFTWrapper<ValueType>::execute( plan, x, x, CUFFT_FORWARD );
    }
    else
    {
        CUFFTWrapper<ValueType>::execute( plan, x, x, CUFFT_INVERSE );
    }
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
    KernelRegistry::set<FFTKernelTrait::fft<RealType<ValueType>>>( CUFFT::fft, ctx, flag );
}

template<>
void CUFFT::RegistratorV<float>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag )
{
}

template<>
void CUFFT::RegistratorV<double>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag )
{
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

CUFFT::CUFFT()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_TYPELIST( SCAI_NUMERIC_TYPES_CUDA )>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ADD );
}

CUFFT::~CUFFT()
{
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_TYPELIST( SCAI_NUMERIC_TYPES_CUDA )>::registerKernels(
        kregistry::KernelRegistry::KERNEL_ERASE );
}

CUFFT CUFFT::guard;    // guard variable for registration

} /* end namespace utilskernel */

} /* end namespace scai */
