/**
 * @file FFTW3Wrapper.hpp
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
 * @brief Wrapper for FFTW3 functions
 * @author Eric Schricker
 * @date 29.09.2016
 */

#pragma once

// local library

// internal scai libraries
#include <scai/common/macros/unused.hpp>
#include <scai/common/macros/assert.hpp>
#include <scai/common/SCAITypes.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/config.hpp>
#include <scai/common/Complex.hpp>

#include <fftw3.h>

#define FFTW3_NAME( name, prefix ) prefix##_##name

namespace scai
{

namespace utilskernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT FFTW3Wrapper;

#define FFTW3WRAPPER_DEF( ValueType, prefix, ComplexType )                                                      \
    template<>                                                                                                  \
    class COMMON_DLL_IMPORTEXPORT FFTW3Wrapper<ValueType>                                                       \
    {                                                                                                           \
    public:                                                                                                     \
        typedef int FFTW3IndexType;                                                                             \
        typedef unsigned int FFTW3FlagType;                                                                     \
        typedef FFTW3_NAME( plan, prefix ) FFTW3PlanType;                                                       \
                                                                                                                \
        static FFTW3PlanType plan_dft_1d( const FFTW3IndexType n,                                               \
            const common::Complex<ValueType> in[],                                                              \
            common::Complex<ValueType> out[],                                                                   \
            const FFTW3IndexType direction,                                                                     \
            const FFTW3FlagType type )                                                                          \
        {                                                                                                       \
            typedef common::Complex<ValueType> SCAIComplex;                                                     \
            SCAIComplex* n_in = const_cast<SCAIComplex*>( in );                                                 \
            ComplexType* in_ = reinterpret_cast<ComplexType*>( n_in );                                          \
            ComplexType* out_ = reinterpret_cast<ComplexType*>( out );                                          \
            return FFTW3_NAME( plan_dft_1d, prefix )( n, in_, out_, direction, type );                          \
        }                                                                                                       \
                                                                                                                \
        static void execute( FFTW3PlanType p )                                                                  \
        {                                                                                                       \
            FFTW3_NAME( execute, prefix )( p );                                                                 \
        }                                                                                                       \
                                                                                                                \
        static void destroy_plan( FFTW3PlanType p )                                                             \
        {                                                                                                       \
            FFTW3_NAME( destroy_plan, prefix )( p );                                                            \
        }                                                                                                       \
    };


FFTW3WRAPPER_DEF( float, fftwf, fftwf_complex );
FFTW3WRAPPER_DEF( double, fftw, fftw_complex );
FFTW3WRAPPER_DEF( long double, fftwl, fftwl_complex );

#ifdef SCAI_COMPLEX_SUPPORTED
//FFTW3WRAPPER_DEF( ComplexFloat, c, sc, ic, dotc );
//FFTW3WRAPPER_DEF( ComplexDouble, z, dz, iz, dotc );
#endif

#undef FFTW3WRAPPER_DEF

} /* end namespace utilskernel */

} /* end namespace scai */

