/**
 * @file Walltime.cpp
 *
 * @license
 * Copyright (c) 2009-2015
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
 * @date 11.06.2015
 */

#include <tracing/Walltime.hpp>

#if defined( _OPENMP )

// with OpenMP support the routine omp_get_wtime can be used

#include <omp.h>

#endif

#if defined( _WIN32 )

#include <windows.h>

#else

#include <sys/time.h>

#endif

namespace common
{

INTEGER_8 Walltime::timestamp()
{
#if defined( WIN32 )

    SYSTEMTIME lpSystemTime;
    GetLocalTime( &lpSystemTime );

    INTEGER_8 ticks = lpSystemTime.wHour;
    ticks = ticks * 60 + lpSystemTime.wMinute;
    ticks = ticks * 60 + lpSystemTime.wSecond;
    ticks = ticks * 1000 + lpSystemTime.wMilliseconds;

#else

    struct timeval tp;
    struct timezone tzp;

    gettimeofday( &tp, &tzp );

    INTEGER_8 ticks = tp.tv_sec;
    ticks = ticks * 1000000 + tp.tv_usec;

#endif

    return ticks;
}

INTEGER_8 Walltime::timerate()
{
#if defined( WIN32 )
    return static_cast<INTEGER_8>( 1000 );       // one tick is a ms
#else
    return static_cast<INTEGER_8>( 1000000 );    // one tick is a us
#endif
}

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

