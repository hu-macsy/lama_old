/**
 * @file utilskernel/cuda/CUFFTWrapper.hpp
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
 * @brief ToDo: Missing description in ./blaskernel/cuda/CUBLASWrapper.hpp
 * @author eschricker
 * @date 24.08.2015
 */

#pragma once

// internal scai libraries
#include <scai/common/cuda/CUDAError.hpp>
#include <scai/common/TypeTraits.hpp>
#include <scai/common/Complex.hpp>

// CUDA
#ifdef SCAI_COMPLEX_SUPPORTED
#include <cuComplex.h>
#endif

#include <cufft.h>

#define CUFFT_NAME( name, P1, P2 ) cufft##name##P1##2##P2

namespace scai
{

namespace utilskernel
{

template<typename ValueType>
class COMMON_DLL_IMPORTEXPORT CUFFTWrapper;

#define CUFFTWRAPPER_DEF( ValueType, ComplexType, R, C )                                            \
    template<>                                                                                      \
    class COMMON_DLL_IMPORTEXPORT CUFFTWrapper<ValueType>                                           \
    {                                                                                               \
    public:                                                                                         \
                                                                                                    \
        typedef int CUFFTDirection;                                                                 \
                                                                                                    \
        static cufftType getTypeC2C()                                                               \
        {                                                                                           \
            return CUFFT_##C##2##C;                                                                 \
        }                                                                                           \
                                                                                                    \
        static cufftType getTypeR2C()                                                               \
        {                                                                                           \
            return CUFFT_##R##2##C;                                                                 \
        }                                                                                           \
                                                                                                    \
        static cufftType getTypeC2R()                                                               \
        {                                                                                           \
            return CUFFT_##C##2##R;                                                                 \
        }                                                                                           \
                                                                                                    \
        /* C2C */                                                                                   \
        static void execute(                                                                        \
                    cufftHandle plan,                                                               \
                    const common::Complex<ValueType> in[],                                          \
                    common::Complex<ValueType> out[],                                               \
                    CUFFTDirection direction )                                                      \
        {                                                                                           \
            const ComplexType* _in = reinterpret_cast<const ComplexType*>( in );                    \
            ComplexType* __in = const_cast<ComplexType*>(_in);                                      \
            ComplexType* _out = reinterpret_cast<ComplexType*>( out );                              \
            SCAI_CUFFT_CALL( CUFFT_NAME( Exec, C, C )( plan, __in, _out, direction ), "execute" )   \
        }                                                                                           \
                                                                                                    \
        /* R2C */                                                                                   \
        static void execute(                                                                        \
                    cufftHandle plan,                                                               \
                    const ValueType in[],                                                           \
                    common::Complex<ValueType> out[] )                                              \
        {                                                                                           \
            ComplexType* _out = reinterpret_cast<ComplexType*>( out );                              \
            ValueType* _in = const_cast<ValueType*>( in );                                          \
            SCAI_CUFFT_CALL( CUFFT_NAME( Exec, R, C )( plan, _in, _out ), "execute ")               \
        }                                                                                           \
                                                                                                    \
        /* C2R */                                                                                   \
        static void execute(                                                                        \
                    cufftHandle plan,                                                               \
                    const common::Complex<ValueType> in[],                                          \
                    ValueType out[] )                                                               \
        {                                                                                           \
            const ComplexType* _in = reinterpret_cast<const ComplexType*>( in );                    \
            ComplexType* __in = const_cast<ComplexType*>( _in );                                    \
            SCAI_CUFFT_CALL( CUFFT_NAME( Exec, C, R )( plan, __in, out ), "execute" )               \
        }                                                                                           \
                                                                                                    \
    };


#ifdef SCAI_COMPLEX_SUPPORTED
    CUFFTWRAPPER_DEF( float, cufftComplex, R, C )
    CUFFTWRAPPER_DEF( double, cufftDoubleComplex, D, Z )
#endif

#undef CUFFT_NAME
#undef CUFFTWRAPPER_DEF

} /* end namespace utilskernel */

} /* end namespace scai */

