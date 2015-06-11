/**
 * @file TraceRegionRecord.cpp
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
 * @brief Implementation of TraceRegionRecord class.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <tracing/TraceRegionRecord.hpp>

#include <tracing/Walltime.hpp>

// others
#include <tracing/RegionTable.hpp>
#include <tracing/TraceConfig.hpp>
#include <tracing/CallTree.hpp>
#include <tracing/VTInterface.hpp>

#include <cstdio>

namespace tracing
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( TraceRegionRecord::logger, "TraceRegionRecord" )

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::start( const char* regionName, const char* file, int lno )
{
    // enter a region withtout using member variables

    if( !TraceConfig::globalTraceFlag )
    {
        return;
    }

    boost::shared_ptr<TraceConfig> traceConfig = TraceConfig::getInstancePtr();

    if( !traceConfig->isEnabled() )
    {
        return;
    }

    double startTime = lama::Walltime::get();

    RegionTable* regionTable = traceConfig->getRegionTable();

    if( !regionTable )
    {
        return; // Tracing is disabled here
    }

    int regionId = regionTable->getRegion( regionName, file, lno );

    LAMA_LOG_DEBUG( logger,
                    "Thread " << regionTable->getId() << ": enters region " << regionName << ", timer = " << regionId )

    if( traceConfig->isTimeTraceEnabled() )
    {
        regionTable->start( regionId, startTime );
    }

    if( traceConfig->isCallTreeEnabled() )
    {
        CallTree::enter( regionId, regionTable->getRegion( regionId ), startTime );
    }

    if( traceConfig->isVampirTraceEnabled() )
    {
        VTInterface::enter( regionTable->getRegion( regionId ) );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::stop( const char* regionName )
{
    // enter a region withtout using member variables

    if( !TraceConfig::globalTraceFlag )
    {
        return;
    }

    boost::shared_ptr<TraceConfig> traceConfig = TraceConfig::getInstancePtr();

    if( !traceConfig->isEnabled() )
    {
        return;
    }

    RegionTable* regionTable = traceConfig->getRegionTable();

    if( !regionTable )
    {
        return; // Tracing is disabled here
    }

    // getCurrentRegionId very fast, matches name against call stack if available

    int regionId = regionTable->getCurrentRegionId( regionName );

    double stopTime = lama::Walltime::get();

    LAMA_LOG_DEBUG( logger, "Thread " << regionTable->getId() << ": leaves region " << regionName )

    if( traceConfig->isTimeTraceEnabled() )
    {
        regionTable->stop( regionId, stopTime );
    }

    if( traceConfig->isCallTreeEnabled() )
    {
        CallTree::leave( regionId, regionTable->getRegion( regionId ), stopTime );
    }

    if( traceConfig->isVampirTraceEnabled() )
    {
        VTInterface::leave( regionTable->getRegion( regionId ) );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::enter( const char* regionName, const char* file, int lno )
{
    mStartTime = lama::Walltime::get();

    mRegionTable = mTraceConfig->getRegionTable();

    if( !mRegionTable )
    {
        return; // Tracing is disabled here
    }

    mRegionId = mRegionTable->getRegion( regionName, file, lno );

    LAMA_LOG_DEBUG( logger,
                    "Thread " << mRegionTable->getId() << ": enters region " << regionName << ", timer = " << mRegionId )

    mTimeTrace = mTraceConfig->isTimeTraceEnabled();

    if( mTimeTrace )
    {
        mRegionTable->start( mRegionId, mStartTime );
    }

    mCallTree = mTraceConfig->isCallTreeEnabled();

    if( mCallTree )
    {
        CallTree::enter( mRegionId, mRegionTable->getRegion( mRegionId ), mStartTime );
    }

    mVampirTrace = mTraceConfig->isVampirTraceEnabled();

    if( mVampirTrace )
    {
        VTInterface::enter( mRegionTable->getRegion( mRegionId ) );
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, const char* file, int lno )
{
    mRegionTable = NULL;

    if( !TraceConfig::globalTraceFlag )
    {
        // tracing is switched off in source code
        return;
    }

    mTraceConfig = TraceConfig::getInstancePtr();

    if( !mTraceConfig->isEnabled() )
    {
        return;
    }

    enter( regionName, file, lno );
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, int n, const char* file, int lno )
{
    mRegionTable = NULL;

    if( !TraceConfig::globalTraceFlag )
    {
        // tracing is switched off in source code
        return;
    }

    mTraceConfig = TraceConfig::getInstancePtr();

    if( !mTraceConfig->isEnabled() )
    {
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
    if( !mRegionTable )
    {
        // tracing was not enabled at all

        return;
    }

    double stopTime = lama::Walltime::get();

    LAMA_LOG_DEBUG( logger,
                    "Thread " << mRegionTable->getId() << ": leaves region " << mRegionTable->getRegion( mRegionId ).getRegionName() << ", timer = " << mRegionId )

    if( mTimeTrace )
    {
        mRegionTable->stop( mRegionId, stopTime );
    }

    if( mCallTree )
    {
        CallTree::leave( mRegionId, mRegionTable->getRegion( mRegionId ), stopTime );
    }

    if( mVampirTrace )
    {
        VTInterface::leave( mRegionTable->getRegion( mRegionId ) );
    }
}

/* -------------------------------------------------------------------------- */

double TraceRegionRecord::spentLast( const char* name )
{
    RegionTable* regionTable = TraceConfig::getInstance().getRegionTable();

    if( !regionTable )
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
