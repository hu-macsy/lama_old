/**
 * @file OpenMP.hpp
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
 * @brief Common defintions for optional use of OpenMP
 * @author Thomas Brandes
 * @date 11.06.2013
 */

#pragma once

#ifdef _OPENMP

#include <omp.h>

#else

#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_set_num_threads( x )

#if defined( WIN32 )

#include <windows.h>

inline double omp_get_wtime( void )
{
    SYSTEMTIME lpSystemTime;
    GetLocalTime( &lpSystemTime );
    return ( lpSystemTime.wHour * 60.0 + lpSystemTime.wMinute ) * 60.0 +
           lpSystemTime.wSecond + lpSystemTime.wMilliseconds * 0.001;
}

#else

#include <sys/time.h>

inline double omp_get_wtime( void )
{
    struct timeval tp;
    struct timezone tzp;
    gettimeofday( &tp, &tzp );
    return ( double ) tp.tv_sec + tp.tv_usec * 0.000001;
}

#endif

#endif

/** atomicAdd used for reductions as reduction directive is unsupported for complex numbers.
 *
 *  Note: template specialization used for float and double
 */

template<typename ValueType>
inline void atomicAdd( ValueType& sharedResult, const ValueType& threadResult )
{
    #pragma omp critical
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( float& sharedResult, const float& threadResult )
{
    #pragma omp atomic
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( double& sharedResult, const double& threadResult )
{
    #pragma omp atomic
    sharedResult += threadResult;
}

template<>
inline void atomicAdd( int& sharedResult, const int& threadResult )
{
    #pragma omp atomic
    sharedResult += threadResult;
}

/* --------------------------------------------------------------------------- */

template<typename ValueType>
inline ValueType atomicInc( ValueType& var )
{
    ValueType val;

    #pragma omp critical
    {
        val = var++;
    }

    return val;
}

template<>
inline int atomicInc( int& var )
{
    return __sync_fetch_and_add( &var, 1 );
}

/** This routine computes within a parallel region a contiguous block for each thread
 *
 *  @param[in] n the range 0 <= i < n is distributed
 *  @param[out] lb, ub is the range belonging to this thread, lb <= i < ub
 *
 *  This routine corresponds the static worksharing of a for( i = 0; i < N; ++i ) loop
 */

template<typename IndexType>
inline void omp_get_my_range( IndexType& lb, IndexType& ub, const IndexType n )
{
    int rank = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    IndexType blockSize = ( n + nthreads - 1 ) / nthreads;
    lb = rank * blockSize;
    ub = ( rank + 1 ) * blockSize;

    if ( ub >= n )
    {
        ub = n;
    }
}

