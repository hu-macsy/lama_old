/**
 * @file Walltime.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Implementation of method that delivers the walltime
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/common/Walltime.hpp>

#if defined( _WIN32 )

#include <windows.h>

#else

#include <sys/time.h>

#endif

namespace scai
{

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
#if defined( WIN32 )

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

} /* end namespace common */

} /* end namespace scai */
