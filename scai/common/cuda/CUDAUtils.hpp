/**
 * @file common/cuda/CUDAUtils.hpp
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
 * @brief Definition of struct that contains util routines overloaded for all possible index types
 * @author Thomas Brandes
 * @date 15.09.2016
 */

#pragma once

#include <cuda.h>

extern "C"
{
    extern __device__ __device_builtin__ int          __iAtomicAdd( int* address, int val );
    extern __device__ __device_builtin__ unsigned int __uAtomicAdd( unsigned int* address, unsigned int val );
    extern __device__ __device_builtin__ unsigned long long __ullAtomicAdd( unsigned long long* address, unsigned long long val );
    extern __device__ __device_builtin__ int          __iAtomicCAS( int* address, int compare, int val );
    extern __device__ __device_builtin__ unsigned int __uAtomicCAS( unsigned int* address, unsigned int compare, unsigned int val );
    extern __device__ __device_builtin__ unsigned long long __ullAtomicCAS( unsigned long long* address, unsigned long long compare, unsigned long long val );
}

namespace scai
{

namespace common
{

/** Structure (instead of namespace) that contains the required utility functions
 *  for each possible index type.
 */

struct CUDAUtils
{
    // atomic compare and swap

    static __inline__ __device__ int atomicCAS( int* address, int compare, int val );
    static __inline__ __device__ unsigned int atomicCAS( unsigned int* address, unsigned int compare, unsigned int val );
    static __inline__ __device__ long atomicCAS( long* address, long compare, long val );
    static __inline__ __device__ unsigned long atomicCAS( unsigned long* address, unsigned long compare, unsigned long val );

    // atomic add

    static __inline__ __device__ int atomicAdd( int* address, int val );
    static __inline__ __device__ unsigned int atomicAdd( unsigned int* address, unsigned int val );
    static __inline__ __device__ long atomicAdd( long* address, long val );
    static __inline__ __device__ unsigned long atomicAdd( unsigned long* address, unsigned long val );

    static __inline__ __device__ void atomicAdd( float* address, float val );
    static __inline__ __device__ void atomicAdd( double* address, double val );

#ifdef SCAI_COMPLEX_SUPPORTED
    static __inline__ __device__ void atomicAdd( ComplexFloat* address, ComplexFloat val );
    static __inline__ __device__ void atomicAdd( ComplexDouble* address, ComplexDouble val );
#endif

};

// -------------------------------- atomicCAS --------------------------------------

__device__ int CUDAUtils::atomicCAS( int* address, int compare, int val )
{
    return __iAtomicCAS( address, compare, val );
}

__device__ unsigned int CUDAUtils::atomicCAS( unsigned int* address, unsigned int compare, unsigned int val )
{
    return __uAtomicCAS( address, compare, val );
}

__device__ long CUDAUtils::atomicCAS( long* address, long compare, long val )
{
    typedef unsigned long long int RepT;

    RepT* ptrCompare = reinterpret_cast<RepT*>( &compare );
    RepT* ptrVal     = reinterpret_cast<RepT*>( &val );

    RepT* t_address = reinterpret_cast<RepT*>( address );

    return __ullAtomicCAS( t_address, *ptrCompare, *ptrVal );
}

__device__ unsigned long CUDAUtils::atomicCAS( unsigned long* address, unsigned long compare, unsigned long val )
{
    typedef unsigned long long int RepT;

    RepT* ptrCompare = reinterpret_cast<RepT*>( &compare );
    RepT* ptrVal     = reinterpret_cast<RepT*>( &val );

    RepT* t_address = reinterpret_cast<RepT*>( address );

    return __ullAtomicCAS( t_address, *ptrCompare, *ptrVal );
}

// -------------------------------- atomicAdd --------------------------------------

__device__ int CUDAUtils::atomicAdd( int* address, int val )
{
    return __iAtomicAdd( address, val );
}

__device__ unsigned int CUDAUtils::atomicAdd( unsigned int* address, unsigned int val )
{
    return __uAtomicAdd( address, val );
}

__device__ unsigned long CUDAUtils::atomicAdd( unsigned long* address, unsigned long val )
{
    typedef unsigned long long int RepT;

    RepT* ptrVal    = reinterpret_cast<RepT*>( &val );
    RepT* t_address = reinterpret_cast<RepT*>( address );

    return __ullAtomicAdd( t_address, *ptrVal );
}

__device__ long CUDAUtils::atomicAdd( long* address, long val )
{
    typedef unsigned long long int RepT;

    RepT* ptrVal    = reinterpret_cast<RepT*>( &val );
    RepT* t_address = reinterpret_cast<RepT*>( address );

    return __ullAtomicAdd( t_address, *ptrVal );
}

__device__ void CUDAUtils::atomicAdd( double* address, double val )
{
    unsigned long long int* address_as_ull =
        ( unsigned long long int* )address;
    unsigned long long int old = *address_as_ull, assumed;

    do
    {
        assumed = old;
        old = __ullAtomicCAS( address_as_ull,
                              assumed,
                              __double_as_longlong( val + __longlong_as_double( assumed ) ) );

        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)

    }
    while ( assumed != old );
}

__device__ void CUDAUtils::atomicAdd( float* address, float val )

{
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 200
    // CUDA runtime offers faster solution for capability >= 2.0
    ::atomicAdd( address, val );
#else
    // old slow solution
    int i_val = __float_as_int( val );
    int tmp0 = 0;
    int tmp1;

    while ( ( tmp1 = atomicCAS( ( int* ) address, tmp0, i_val ) ) != tmp0 )
    {
        tmp0 = tmp1;
        i_val = __float_as_int( val + __int_as_float( tmp1 ) );
    }

#endif
}

#ifdef SCAI_COMPLEX_SUPPORTED

__device__ void CUDAUtils::atomicAdd( ComplexFloat* address, ComplexFloat val )
{
    float* faddress = ( float* ) address;
    atomicAdd( &faddress[0], val.real() );
    atomicAdd( &faddress[1], val.imag() );
}

__device__ void CUDAUtils::atomicAdd( ComplexDouble* address, ComplexDouble val )
{
    double* daddress = ( double* ) address;
    atomicAdd( &daddress[0], val.real() );
    atomicAdd( &daddress[1], val.imag() );
}

#endif

} /* end namespace common */

} /* end namespace scai */
