/**
 * @file Walltime.cpp
 *
 * @license
 * Copyright (c) 2009-2013
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Implementation of method that delivers the walltime
 * @author Thomas Brandes
 * @date 25.04.2013
 * $Id$
 */

// hpp
#include <lama/Walltime.hpp>

#if defined( _OPENMP )

// with OpenMP support the routine omp_get_wtime can be used

#include <omp.h>

#elif defined( _WIN32 )

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

