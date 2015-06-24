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

// others
#include <tracing/RegionTable.hpp>
#include <tracing/TraceConfig.hpp>
#include <tracing/CallTree.hpp>
#include <tracing/VTInterface.hpp>

#include <common/Walltime.hpp>

#include <cstdio>

namespace tracing
{

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( TraceRegionRecord::logger, "TraceRegionRecord" )

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::initSettings()
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

    mRegionTable = mTraceConfig->getRegionTable();

    if ( mRegionTable )
    {
        // get detailed info about what to trace

        mTimeTrace   = mTraceConfig->isTimeTraceEnabled();
        mCallTree    = mTraceConfig->isCallTreeEnabled();
        mVampirTrace = mTraceConfig->isVampirTraceEnabled();
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, const char* file, int lno )
{
    initSettings();

    // if tracing is enable we get the region, may be it will be defined

    if ( mRegionTable )
    {
        mRegionId = mRegionTable->getRegion( regionName, file, lno );
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName )
{
    initSettings();

    if ( mRegionTable )
    {
        // getCurrentRegionId very fast, matches name against call stack if available

        mRegionId = mRegionTable->getCurrentRegionId( regionName );
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, int n, const char* file, int lno )
{
    initSettings();

    // if tracing is enable we get the region, may be it will be defined

    if ( mRegionTable )
    {
        // compose a new name for the region with n as suffix

        std::ostringstream fullRegionName;
    
        fullRegionName << regionName << "_" << n;

        mRegionId = mRegionTable->getRegion( fullRegionName.str().c_str(), file, lno );

        // full region name is no longer needed, will be available by region table
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::~TraceRegionRecord()
{
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::enter()
{
    if( !mRegionTable )
    {
        return;   // tracing disabled 
    }

    CounterArray startCounterValues( true );  // get stamp of all counters

    if( mTimeTrace )
    {
        mRegionTable->start( mRegionId, startCounterValues );
    }

    if( mCallTree )
    {
        CallTree::enter( mRegionId, mRegionTable->getRegion( mRegionId ), startCounterValues );
    }

    if( mVampirTrace )
    {
        VTInterface::enter( mRegionTable->getRegion( mRegionId ) );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::leave()
{
    if( !mRegionTable )
    {
        return;   // tracing was not enabled at all
    }

    CounterArray stopCounterValues( true ); // stamp it

    LAMA_LOG_DEBUG( logger,
                    "Thread " << mRegionTable->getId() << ": leaves region " << mRegionTable->getRegion( mRegionId ).getRegionName()
                     << ", counters = " << stopCounterValues )

    if( mTimeTrace )
    {
        mRegionTable->stop( mRegionId, stopCounterValues );
    }

    if( mCallTree )
    {
        CallTree::leave( mRegionId, mRegionTable->getRegion( mRegionId ), stopCounterValues );
    }

    if( mVampirTrace )
    {
        VTInterface::leave( mRegionTable->getRegion( mRegionId ) );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::start( const char* regionName, const char* file, int lno )
{
    // static version has to build a new record 

    TraceRegionRecord record( regionName, file, lno );
    record.enter();
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::stop( const char* regionName )
{
    // static version has to build a new record 

    TraceRegionRecord record( regionName );  // region should be known here
    record.leave();
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
