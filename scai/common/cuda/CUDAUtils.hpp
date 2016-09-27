/**
 * @file CUDAUtils.hpp
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
 * @brief Definition of struct that contains util routines overloaded for all possible index types
 * @author Thomas Brandes
 * @date 15.09.2016
 */

#pragma once

#include <cuda.h>

namespace scai
{

namespace common
{

/** Structure (instead of namespace) that contains the required utility functions
 *  for each possible index type.
 */

struct CUDAUtils
{
    static __inline__ __device__ int atomicCAS( int* address, int compare, int val );
    static __inline__ __device__ unsigned int atomicCAS( unsigned int* address, unsigned int compare, unsigned int val );
    static __inline__ __device__ long atomicCAS( long* address, long compare, long val );
    static __inline__ __device__ unsigned long atomicCAS( unsigned long* address, unsigned long compare, unsigned long val );

    static __inline__ __device__ int atomicAdd( int* address, int val );
    static __inline__ __device__ unsigned int atomicAdd( unsigned int* address, unsigned int val );
    static __inline__ __device__ long atomicAdd( long* address, long val );
    static __inline__ __device__ unsigned long atomicAdd( unsigned long* address, unsigned long val );
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

    return __ullAtomicCAS( t_address, *ptrCompare, *ptrVal);
}

__device__ unsigned long CUDAUtils::atomicCAS( unsigned long* address, unsigned long compare, unsigned long val )
{
    typedef unsigned long long int RepT;

    RepT* ptrCompare = reinterpret_cast<RepT*>( &compare );
    RepT* ptrVal     = reinterpret_cast<RepT*>( &val );

    RepT* t_address = reinterpret_cast<RepT*>( address );

    return __ullAtomicCAS( t_address, *ptrCompare, *ptrVal);
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

} /* end namespace common */

} /* end namespace scai */
