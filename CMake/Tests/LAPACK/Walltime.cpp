/**
 * @file Walltime.cpp
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
 * @brief Implementation of method that delivers the walltime
 * @author Thomas Brandes
 * @date 25.04.2013
 */

// hpp
#include "Walltime.hpp"

#if defined( _OPENMP )

// with OpenMP support the routine omp_get_wtime can be used

#include <omp.h>

#elif defined( WIN32 )

#include <windows.h>

#else

#include <sys/time.h>  

#endif

namespace lama
{

double Walltime::get()
{
    
#if defined( _OPENMP )

    return omp_get_wtime();

#elif defined( WIN32 )

    SYSTEMTIME lpSystemTime;
    GetLocalTime( &lpSystemTime );
    return ( lpSystemTime.wHour * 60.0 + lpSystemTime.wMinute ) * 60.0 +
           lpSystemTime.wSecond + lpSystemTime.wMilliseconds * 0.001;

#else

    struct timeval tp;
    struct timezone tzp;

    gettimeofday( &tp, &tzp );

    return (double) tp.tv_sec + tp.tv_usec * 0.000001;

#endif 
     
}

} // namespace lama
