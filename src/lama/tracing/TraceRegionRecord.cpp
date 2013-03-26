/**
 * @file TraceRegionRecord.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Implementation of TraceRegionRecord class.
 * @author Thomas Brandes
 * @date 01.09.2011
 * $Id$
 */

// hpp
#include <lama/tracing/TraceRegionRecord.hpp>

// others
#include <lama/tracing/RegionTable.hpp>
#include <lama/tracing/TraceConfig.hpp>
#include <lama/tracing/VTInterface.hpp>

#include <omp.h>

#include <cstdio>

namespace tracing
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( TraceRegionRecord::logger, "TraceRegionRecord" )

/* -------------------------------------------------------------------------- */

static double getWallTime()
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

#endif //WIN32
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::start( const char* regionName, const char* file, int lno )
{
    // enter a region withtout using member variables

    boost::shared_ptr<TraceConfig> traceConfig = TraceConfig::getInstancePtr();

    if ( !traceConfig->isEnabled() )
    {
        return;
    }

    double startTime = getWallTime();

    RegionTable* regionTable = traceConfig->getRegionTable();

    if ( !regionTable )
    {
        return; // Tracing is disabled here
    }

    int regionId = regionTable->getRegion( regionName, file, lno );

    LAMA_LOG_DEBUG( logger,
                    "Thread " << regionTable->getId() << ": enters region " << regionName << ", timer = " << regionId )

    if ( traceConfig->isTimeTraceEnabled() )
    {
        regionTable->start( regionId, startTime );
    }

    if ( traceConfig->isVampirTraceEnabled() )
    {
        VTInterface::enter( regionTable->getRegion( regionId ) );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::stop( const char* regionName )
{
    // enter a region withtout using member variables

    boost::shared_ptr<TraceConfig> traceConfig = TraceConfig::getInstancePtr();

    if ( !traceConfig->isEnabled() )
    {
        return;
    }

    RegionTable* regionTable = traceConfig->getRegionTable();

    if ( !regionTable )
    {
        return; // Tracing is disabled here
    }

    int regionId = regionTable->getCurrentRegionId( regionName );

    double stopTime = getWallTime();

    LAMA_LOG_DEBUG( logger, "Thread " << regionTable->getId() << ": leaves region " << regionName )

    if ( traceConfig->isTimeTraceEnabled() )
    {
        regionTable->stop( regionId, stopTime );
    }

    if ( traceConfig->isVampirTraceEnabled() )
    {
        VTInterface::leave( regionTable->getRegion( regionId ) );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::enter( const char* regionName, const char* file, int lno )
{
    mStartTime = getWallTime();

    mRegionTable = mTraceConfig->getRegionTable();

    if ( !mRegionTable )
    {
        return; // Tracing is disabled here
    }

    mRegionId = mRegionTable->getRegion( regionName, file, lno );

    LAMA_LOG_DEBUG( logger,
                    "Thread " << mRegionTable->getId() << ": enters region " << regionName << ", timer = " << mRegionId )

    mTimeTrace = mTraceConfig->isTimeTraceEnabled();

    if ( mTimeTrace )
    {
        mRegionTable->start( mRegionId, mStartTime );
    }

    mVampirTrace = mTraceConfig->isVampirTraceEnabled();

    if ( mVampirTrace )
    {
        VTInterface::enter( mRegionTable->getRegion( mRegionId ) );
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, const char* file, int lno )
{
    mTraceConfig = TraceConfig::getInstancePtr();

    if ( !mTraceConfig->isEnabled() )
    {
        mRegionTable = NULL;
        return;
    }

    enter( regionName, file, lno );
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, int n, const char* file, int lno )

{
    mTraceConfig = TraceConfig::getInstancePtr();

    if ( !mTraceConfig->isEnabled() )
    {
        mRegionTable = NULL;

        return;
    }

    // compose a new name for the region with n as suffix

    std::ostringstream fullRegionName;

    fullRegionName << regionName << "_" << n;

    enter( fullRegionName.str().c_str(), file, lno );
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::~TraceRegionRecord()
{
    if ( !mRegionTable )
    {
        // tracing was not enabled at all

        return;
    }

    double stopTime = getWallTime();

    LAMA_LOG_DEBUG( logger,
                    "Thread " << mRegionTable->getId() << ": leaves region " << mRegionTable->getRegion( mRegionId ).getRegionName() << ", timer = " << mRegionId )

    if ( mTimeTrace )
    {
        mRegionTable->stop( mRegionId, stopTime );
    }

    if ( mVampirTrace )
    {
        VTInterface::leave( mRegionTable->getRegion( mRegionId ) );
    }
}

/* -------------------------------------------------------------------------- */

double TraceRegionRecord::spentLast( const char* name )
{
    RegionTable* regionTable = TraceConfig::getInstance().getRegionTable();

    if ( !regionTable )
    {
        LAMA_LOG_WARN( logger, "spentLast " << name << ": trace is disabled" )
    }

    int regionId = regionTable->getRegion( name, NULL, 0 );

    const RegionEntry& region = regionTable->getRegion( regionId );

    double lastTime = region.getLastTime();

    LAMA_LOG_DEBUG( logger, "Spent time for last call of " << region.getRegionName() << " : " << lastTime )
    return lastTime;
}

} // namespace
