/**
 * @file TraceRegionRecord.cpp
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
 * @brief Implementation of TraceRegionRecord class.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/TraceRegionRecord.hpp>

// local library
#include <scai/tracing/TraceData.hpp>
#include <scai/tracing/TraceConfig.hpp>
#include <scai/tracing/VTInterface.hpp>

// internal scai libraries
#include <scai/common/Walltime.hpp>

// std
#include <cstdio>

namespace scai
{

namespace tracing
{

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( TraceRegionRecord::logger, "TraceRegionRecord" )

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::initSettings()
{
    mTraceData = NULL;

    if ( !TraceConfig::globalTraceFlag )
    {
        // tracing is switched off in source code
        return;
    }

    mTraceConfig = TraceConfig::getInstancePtr();

    if ( !mTraceConfig->isEnabled() )
    {
        return;
    }

    mTraceData = mTraceConfig->getTraceData();

    if ( mTraceData )
    {
        // get detailed info about what to trace
        mTimeTrace   = mTraceConfig->isTimeTraceEnabled();
        mCallTree    = mTraceConfig->isCallTreeEnabled();
        mVampirTrace = mTraceConfig->isVampirTraceEnabled();
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, const char* file, int scl )
{
    initSettings();

    // if tracing is enable we get the region, may be it will be defined

    if ( mTraceData )
    {
        mRegionId    = mTraceData->getRegionId( regionName, file, scl );
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName )
{
    initSettings();

    if ( mTraceData )
    {
        // getCurrentRegionId very fast, matches name against call stack if available
        mRegionId    = mTraceData->getCurrentRegionId( regionName );
    }
    else
    {
        mRegionId    = 0;
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::TraceRegionRecord( const char* regionName, int n, const char* file, int lno )
{
    initSettings();

    // if tracing is enable we get the region, may be it will be defined

    if ( mTraceData )
    {
        // compose a new name for the region with n as suffix
        std::ostringstream fullRegionName;
        fullRegionName << regionName << "_" << n;
        mRegionId    = mTraceData->getRegionId( fullRegionName.str().c_str(), file, lno );
        // full region name is no longer needed, will be available by region table
    }
    else
    {
        mRegionId    = 0;
    }
}

/* -------------------------------------------------------------------------- */

TraceRegionRecord::~TraceRegionRecord()
{
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::enter()
{
    if ( !mTraceData )
    {
        return;   // tracing disabled
    }

    RegionEntry& regionEntry = mTraceData->getRegion( mRegionId );
    SCAI_LOG_INFO( logger, "enter " << regionEntry.getRegionName() )

    if ( mTimeTrace | mCallTree )
    {
        mTraceData->enter( mRegionId, regionEntry );
    }

    if ( mVampirTrace )
    {
        VTInterface::enter( regionEntry );
    }
}

/* -------------------------------------------------------------------------- */

void TraceRegionRecord::leave()
{
    if ( !mTraceData )
    {
        return;   // tracing was not enabled at all
    }

    RegionEntry& regionEntry = mTraceData->getRegion( mRegionId );
    SCAI_LOG_INFO( logger, "leave " << regionEntry.getRegionName() )

    if ( mTimeTrace | mCallTree )
    {
        mTraceData->leave( mRegionId, regionEntry );
    }

    if ( mVampirTrace )
    {
        VTInterface::leave( regionEntry );
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
    TraceData* traceData = TraceConfig::getInstance().getTraceData();

    if ( !traceData )
    {
        SCAI_LOG_WARN( logger, "spentLast " << name << ": trace is disabled" )
    }

    int regionId = traceData->getRegionId( name, NULL, 0 );
    const RegionEntry& region = traceData->getRegion( regionId );
    double lastTime = region.getLastTime();
    SCAI_LOG_DEBUG( logger, "Spent time for last call of " << region.getRegionName() << " : " << lastTime )
    return lastTime;
}

} /* end namespace tracing */

} /* end namespace scai */
