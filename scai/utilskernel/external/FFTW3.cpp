/**
 * @file FFTW3.cpp
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
 * @brief Implementations that uses FFTW3
 * @author Eric Schricker
 * @date 29.09.2016
 */

// hpp
#include <scai/utilskernel/external/FFTW3.hpp>

// local library
#include <scai/utilskernel/FFTKernelTrait.hpp>
#include <scai/utilskernel/external/FFTW3Types.hpp>
#include <scai/utilskernel/external/FFTW3Wrapper.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/TypeTraits.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Settings.hpp>
#include <scai/common/macros/unused.hpp>


namespace scai
{

using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( FFTW3::logger, "external.FFTW3" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTW3::fft( common::Complex<ValueType> x[], IndexType nb, IndexType n, IndexType m, int dir )
{
    SCAI_REGION( "FFTW3.fft" )

    SCAI_LOG_INFO( logger, "fft<" << common::TypeTraits<ValueType>::id() << "> @ FFTW, " 
                   << nb << " x " << n << " = 2 ** " << m )

    typedef typename FFTW3Wrapper<ValueType>::FFTW3PlanType PlanType;
    typedef typename FFTW3Wrapper<ValueType>::FFTW3IndexType FFTW3IndexType;

    for ( IndexType i = 0; i < nb; ++i )
    {
        PlanType p = NULL;

        p = FFTW3Wrapper<ValueType>::plan_dft_1d( static_cast<FFTW3IndexType>( n ), 
                                                  x + i * n, x + i * n, 
                                                  dir == 1 ? FFTW_FORWARD : FFTW_BACKWARD,
                                                  FFTW_ESTIMATE );

        FFTW3Wrapper<ValueType>::execute( p );

        FFTW3Wrapper<ValueType>::destroy_plan( p );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTW3::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    bool useFFTW = true;
    common::Settings::getEnvironment( useFFTW, "SCAI_USE_FFTW" );

    if ( !useFFTW )
    {
        return;
    }

    const common::ContextType ctx = common::ContextType::Host;
    using kregistry::KernelRegistry;
    SCAI_LOG_INFO( logger,
                   "register FFTW3-routines for Host at kernel registry [" << flag << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<FFTKernelTrait::fft<ValueType> >( fft, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

FFTW3::FFTW3()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_REPLACE;
    kregistry::mepr::RegistratorV<RegistratorV,SCAI_NUMERIC_TYPES_FFTW3_LIST>::registerKernels( flag );
}

FFTW3::~FFTW3()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegistratorV,SCAI_NUMERIC_TYPES_FFTW3_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

FFTW3 FFTW3::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
