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

#include <scai/common/OpenMP.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/Complex.hpp>
#include <scai/common/macros/unused.hpp>


namespace scai
{

using common::TypeTraits;
using common::Complex;

namespace utilskernel
{

SCAI_LOG_DEF_LOGGER( OpenMPFFT::logger, "OpenMP.FFT" )

/* --------------------------------------------------------------------------- */

template<typename ValueType>
static void DFT( int dir, IndexType m, Complex<ValueType> x1[] )
{
    ValueType PI_2 = 2.0 * 3.1415926535;

    std::unique_ptr<Complex<ValueType>[]> x2( new Complex<ValueType>[m] );

    for ( IndexType i = 0; i < m; i++ )
    {
        x2[i] = 0;

        ValueType arg = - dir * PI_2 * ValueType( i ) / ValueType( m );

        for ( IndexType k = 0; k < m; k++ )
        {
            ValueType cosarg = std::cos( k * arg );
            ValueType sinarg = std::sin( k * arg );
            x2[i] += x1[k] * Complex<ValueType>( cosarg , sinarg );
        }
    }

    /* Copy the data back */

    for ( IndexType i = 0; i < m; i++ )
    {
        x1[i] = x2[i];
    }
}

template<typename ValueType>
void OpenMPFFT::fft( Complex<ValueType> x[], IndexType n, IndexType m, int dir )
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

        c.imag( sqrt( ( 1.0 - c.real() ) / 2.0 ) );

        if ( dir == 1 )
        {
            c.imag( -c.imag() );
        }

        c.real( sqrt( ( 1.0 + c.real() ) / 2.0 ) );
    }

    return;
}

template<typename ValueType>
static void FFT( int dir, IndexType m, Complex<ValueType> x[] )
 {
    /*Calculate the number of points */

    IndexType n = 1;

    for ( IndexType i = 0; i < m; i++ )
    {
        n <<= 1;
    }

    OpenMPFFT::fft( x, n, m, dir );
}

bool Powerof2(int n,int *m,int *twopm)
{
   if (n <= 1) {
      *m = 0;
      *twopm = 1;
      return false;
   }

   *m = 1;
   *twopm = 2;
   do {
      (*m)++;
      (*twopm) *= 2;
   } while (2*(*twopm) <= n);

   if (*twopm != n)
      return false;
   else
      return true;
}

template<typename ValueType>
void FFT2D(Complex<ValueType> c[], int nx, int ny, int dir)
{
   int i,j;
   int m,twopm;
   double *real,*imag;

   /* Transform the rows */

    std::unique_ptr<Complex<ValueType>[]> col( new Complex<ValueType>[nx] );

    if (!Powerof2(nx,&m,&twopm) || twopm != nx)
      return;

   for (j=0;j<ny;j++) {
      for (i=0;i<nx;i++) {
         col[i] = c[i * ny + j];
      }
      FFT(dir,m,col);
      for (i=0;i<nx;i++) {
         c[i * ny + j] = col[i];
      }
   }

   /* Transform the columns */

   if (!Powerof2(ny,&m,&twopm) || twopm != ny)
      return;
   for (i=0;i<nx;i++) {
      FFT(dir,m,&c[ i * ny ] );
   }
}

template<typename ValueType>
void OpenMPFFT::paddedForward1D(
    const IndexType n,
    const IndexType npad,
    const ValueType in[],
    Complex<ValueType> out[] )
{
    // compute m, n2 = 2**m with n2 >= n

    IndexType m = 0;
    IndexType n2 = 1;

    while ( n2 < n || n2 < npad )
    {
        n2 *= 2;
        m++;
    }

    SCAI_LOG_ERROR( logger, "FFT1D_forward, n = " << n << ", npad = " << npad << ", n2 = " << n2 << " = 2 ** " << m )

    std::unique_ptr<Complex<ValueType>[]> arr( new Complex<ValueType>[n2] );

    for ( IndexType i = 0; i < n; ++i )
    {
        arr[i] = in[i];
    }

    for ( IndexType i = n; i < n2; ++i )
    {
        arr[i] = 0;
    }

    // DFT( 1, n2, arr.get()
    FFT( 1, m, arr.get() );

    for ( IndexType i = 0; i < npad; ++i )
    {
        out[i] = arr[i];
    }
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
void OpenMPFFT::paddedBackward1D(
    const IndexType n,
    const IndexType npad,
    const ValueType in[],
    common::Complex<ValueType> out[] )
{
    // compute m, n2 = 2**m with n2 >= n

    IndexType m = 0;
    IndexType n2 = 1;

    while ( n2 < n || n2 < npad )
    {
        n2 <<= 1;
        m++;
    }

    SCAI_LOG_ERROR( logger, "FFT1D_backward, n = " << n << ", npad = " << npad << ", n2 = " << n2 << " = 2 ** " << m )

    std::unique_ptr<Complex<ValueType>[]> arr( new Complex<ValueType>[n2] );

    for ( IndexType i = 0; i < n; ++i )
    {
        arr[i] = in[i];
    }

    for ( IndexType i = n; i < n2; ++i )
    {
        arr[i] = 0;
    }

    // DFT( 1, n2, arr.get()
    FFT( -1, m, arr.get() );

    for ( IndexType i = 0; i < npad; ++i )
    {
        out[i] = arr[i];
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
    // KernelRegistry::set<FFTKernelTrait::paddedForward1D<ValueType> >( paddedForward1D, ctx, flag );
    // KernelRegistry::set<FFTKernelTrait::paddedBackward1D<ValueType> >( paddedBackward1D, ctx, flag );
    KernelRegistry::set<FFTKernelTrait::fft<ValueType> >( fft, ctx, flag );
}

/* --------------------------------------------------------------------------- */
/*    Constructor/Desctructor with registration                                */
/* --------------------------------------------------------------------------- */

#define SCAI_NUMERIC_TYPES_FFT float, double

#define SCAI_NUMERIC_TYPES_FFT_LIST SCAI_TYPELIST( SCAI_NUMERIC_TYPES_FFT )

OpenMPFFT::OpenMPFFT()
{
    SCAI_LOG_INFO( logger, "register UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ADD;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_FFT_LIST>::registerKernels( flag );
}

OpenMPFFT::~OpenMPFFT()
{
    SCAI_LOG_INFO( logger, "unregister UtilsKernel OpenMP-routines for Host" )
    const kregistry::KernelRegistry::KernelRegistryFlag flag = kregistry::KernelRegistry::KERNEL_ERASE;
    kregistry::mepr::RegistratorV<RegistratorV, SCAI_NUMERIC_TYPES_FFT_LIST>::registerKernels( flag );
}

/* --------------------------------------------------------------------------- */
/*    Static variable to force registration during static initialization      */
/* --------------------------------------------------------------------------- */

OpenMPFFT OpenMPFFT::guard;

} /* end namespace utilskernel */

} /* end namespace scai */
