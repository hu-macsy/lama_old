/**
 * @file OpenMPFFT.cpp
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
 * @brief Implementations that uses OpenMPFFT
 * @author Eric Schricker
 * @date 29.09.2016
 */

// hpp
#include <scai/utilskernel/openmp/OpenMPFFT.hpp>
#include <scai/utilskernel/FFTKernelTrait.hpp>

// internal scai libraries
#include <scai/kregistry/KernelRegistry.hpp>

#include <scai/tracing.hpp>

#include <scai/common/TypeTraits.hpp>
#include <scai/common/Math.hpp>

#include <scai/common/OpenMP.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/macros/unused.hpp>

namespace scai
{

using common::TypeTraits;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPFFT::logger, "OpenMP.FFT" )

/* --------------------------------------------------------------------------- */

#ifdef SCAI_COMPLEX_SUPPORTED

using common::Complex;

template<typename ValueType>
void OpenMPFFT::fft1( Complex<ValueType> x[], IndexType n, IndexType m, int dir )
{
    IndexType i, i1, i2, j, k, l, l1, l2;
    Complex<ValueType> tx, t1, u, c;

    /* Do the bit reversal */

    i2 = n >> 1;
    j = 0;

    for ( i = 0; i < n - 1 ; i++ )
    {
        if ( i < j )
        {
            std::swap( x[i], x[j] );
        }

        k = i2;

        while ( k <= j )
        {
            j -= k;
            k >>= 1;
        }

        j += k;
    }

    /* Compute the FFT */
    c = Complex<ValueType>( -1, 0 );
    l2 = 1;

    for ( l = 0; l < m; l++ )
    {
        l1 = l2;
        l2 <<= 1;
        u.real( 1.0 );
        u.imag( 0.0 );

        for ( j = 0; j < l1; j++ )
        {
            for ( i = j; i < n; i += l2 )
            {
                i1 = i + l1;
                t1 = u * x[i1];
                x[i1] = x[i] - t1;
                x[i] += t1;
            }

            u = u * c;
        }

        c.imag( common::Math::sqrt( ( 1.0 - c.real() ) / 2.0 ) );

        if ( dir == 1 )
        {
            c.imag( -c.imag() );
        }

        c.real( common::Math::sqrt( ( 1.0 + c.real() ) / 2.0 ) );
    }
}

template<typename ValueType>
void OpenMPFFT::fft( Complex<ValueType> x[], IndexType k, IndexType n, IndexType m, int dir )
{
    SCAI_REGION( "OpenMP.fft" )

    SCAI_ASSERT_EQ_DEBUG( ( 1 << m ), n , " 2 ** m ( " << m << " ) != " << n << " = n" )

    SCAI_LOG_INFO( logger, "fft<" << common::TypeTraits<ValueType>::id() << "> @ OpenMP, "
                     << k << " x " << n << " = 2 ** " << m )

    #pragma omp parallel for
    for ( IndexType i = 0; i < k; ++i )
    {
        OpenMPFFT::fft1( x + i * n, n, m, dir );
    }
}

/* --------------------------------------------------------------------------- */
/*     Template instantiations via registration routine                        */
/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPFFT::RegistratorV<ValueType>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag flag )
{
    const common::ContextType ctx = common::ContextType::Host;
    using kregistry::KernelRegistry;
    SCAI_LOG_INFO( logger,
                   "register OpenMPFFT-routines for Host at kernel registry [" << flag << " --> " << common::getScalarType<ValueType>() << "]" )
    KernelRegistry::set<FFTKernelTrait::fft<RealType<ValueType>> >( fft, ctx, flag );
}

#endif

template<>
void OpenMPFFT::RegistratorV<float>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag )
{
}

template<>
void OpenMPFFT::RegistratorV<double>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag )
{
}

template<>
void OpenMPFFT::RegistratorV<scai::LongDouble>::registerKernels( kregistry::KernelRegistry::KernelRegistryFlag )
{
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

OpenMPFFT::OpenMPFFT()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

OpenMPFFT::~OpenMPFFT()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_HOST_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPFFT OpenMPFFT::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
