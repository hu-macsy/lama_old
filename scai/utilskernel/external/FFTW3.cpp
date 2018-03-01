/**
 * @file FFTW3.cpp
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
#include <scai/common/macros/unused.hpp>


namespace scai
{

using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( FFTW3::logger, "external.FFTW3" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTW3::paddedForward1D(
    const IndexType n,
    const IndexType npad,
    const ValueType in[],
    common::Complex<ValueType> out[] )
{
    SCAI_REGION( "external.FFTW3.paddedForward1D" )

    SCAI_LOG_INFO( logger,
                  "fftw padded forward 1-dimensional with: n(" << n << ") npad(" << npad << ") in(" << in << ") and out(" << out << ")" )

    typedef common::Complex<ValueType> ComplexType;
    typedef typename FFTW3Wrapper<ValueType>::FFTW3PlanType PlanType;
    typedef typename FFTW3Wrapper<ValueType>::FFTW3IndexType FFTW3IndexType;

    std::unique_ptr<ComplexType[]> in1( new ComplexType[npad] );

    PlanType p = NULL;

    p = FFTW3Wrapper<ValueType>::plan_dft_1d( static_cast<FFTW3IndexType>( npad ), in1.get(), out, FFTW_FORWARD,
                                              FFTW_ESTIMATE );

    const ComplexType zero = ComplexType( 0 );

    for( IndexType i = 0; i < n; ++i )
    {
        in1[i] = in[i];
    }

    for( IndexType i = n; i < npad; ++i )
    {
        in1[i] = zero;
    }

    FFTW3Wrapper<ValueType>::execute( p );

    FFTW3Wrapper<ValueType>::destroy_plan( p );
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTW3::paddedBackward1D(
    const IndexType n,
    const IndexType npad,
    const ValueType in[],
    common::Complex<ValueType> out[] )
{
SCAI_REGION( "external.FFTW3.paddedBackward1D" )

        SCAI_LOG_INFO( logger,
                   "fftw padded backward 1-dimensional with: n(" << n << ") npad(" << npad << ") in(" << in << ") and out(" << out << ")" )

    typedef common::Complex<ValueType> ComplexType;
    typedef typename FFTW3Wrapper<ValueType>::FFTW3PlanType PlanType;
    typedef typename FFTW3Wrapper<ValueType>::FFTW3IndexType FFTW3IndexType;

    std::unique_ptr<ComplexType[]> in1( new ComplexType[npad] );

    PlanType p = NULL;

    p = FFTW3Wrapper<ValueType>::plan_dft_1d( static_cast<FFTW3IndexType>( npad ), in1.get(), out, FFTW_BACKWARD,
                                              FFTW_ESTIMATE );

    const ComplexType zero = ComplexType( 0 );

    for( IndexType i = 0; i < n; ++i )
    {
        in1[i] = in[i];
    }

    for( IndexType i = n; i < npad; ++i )
    {
        in1[i] = zero;
    }

    FFTW3Wrapper<ValueType>::execute( p );

    FFTW3Wrapper<ValueType>::destroy_plan( p );
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void FFTW3::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::ContextType ctx = common::ContextType::Host;
    using kregistry::KernelRegistry;
    SCAI_LOG_INFO( logger,
                   "register FFTW3-routines for Host at kernel registry [" << flag << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<FFTKernelTrait::paddedForward1D<ValueType> >( paddedForward1D, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

FFTW3::FFTW3()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
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
