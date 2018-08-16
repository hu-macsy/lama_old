/**
 * @file TraceData.cpp
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
 * @brief Definition of class that contains all data needed for tracing.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

/// local library
#include <scai/tracing/TraceData.hpp>
#include <scai/tracing/CallStack.hpp>
#include <scai/tracing/CallTreeTable.hpp>

// internal scai libraries
#include <scai/common/macros/throw.hpp>
#include <scai/common/macros/assert.hpp>

namespace scai
{

namespace tracing
{

SCAI_LOG_DEF_LOGGER( TraceData::logger, "TraceData" )

void TraceData::enter( const int regionId, RegionEntry& region )
{
//    SCAI_ASSERT( &region != NULL, "NULL pointer for region" )
    CounterArray enterCounterValues( true );  // get stamp of all counters
    SCAI_LOG_DEBUG( logger, "enter " << regionId << ", region= " << &region )
    // SCAI_LOG_DEBUG( logger, "enter " << regionId << ", " << region << ", counters = " << enterCounterValues )

    if ( mCallTreeTable.get() != NULL )
    {
        if ( region.firstAccess() )
        {
            // write Info about the region in the call tree
            mCallTreeTable->writeRegion( regionId, region ) ;
        }

        if ( mCallStack.empty() )
        {
            // This is the first time we enter a region, do initialization
            mCallTreeTable->initCounters( enterCounterValues );
        }
        else
        {
            // Before we enter the called region add exclusive costs so far for caller region
            int caller_region = mCallStack.currentRegionId();
            int scl = 0;
            mCallTreeTable->addExclusiveCosts( caller_region, scl, enterCounterValues );
        }
    }

    mCallStack.push( regionId, enterCounterValues );
};

void TraceData::leave( const int regionId, RegionEntry& region )
{
//    SCAI_ASSERT( &region != NULL, "NULL pointer for region" )
    CounterArray leaveCounterValues( true );  // get stamp of all counters
    SCAI_LOG_DEBUG( logger, "leave " << regionId << ", region = " << &region )

    // SCAI_LOG_DEBUG( logger, "leave " << regionId << ", " << region << ", counters = " << leaveCounterValues )

    if ( mCallStack.empty() )
    {
        SCAI_LOG_ERROR( logger, "stop region on empty call region stack" )
        return;
    }

    const int currentRegionId = mCallStack.currentRegionId();

    SCAI_ASSERT_EQUAL( currentRegionId, regionId,
                       "mismatch call stack, current region = "
                       << mRegionTable.getRegion( currentRegionId ).getRegionName()
                       << ", stop for " << region.getRegionName() )
    double spentTime = leaveCounterValues.getWalltime( mCallStack.currentCounters() );

    CounterArray costs;

    mCallStack.getCosts( costs, leaveCounterValues );  // costs = counterVals - startVals

    region.addCall( spentTime );

    SCAI_LOG_DEBUG( logger, "Region " << regionId << ": spent time = " << spentTime )

    if ( mCallTreeTable.get() != NULL )
    {
        mCallTreeTable->addExclusiveCosts( regionId, 0, leaveCounterValues );
    }

    SCAI_LOG_DEBUG( logger, region.getRegionName() << ", spent time = " << spentTime << ", costs = " << costs )
    mCallStack.pop();

    // correct exclusive time of previous entry in call stack
    if ( !mCallStack.empty() )
    {
        int callRegionId = mCallStack.currentRegionId();
        RegionEntry& callRegion = mRegionTable.getRegion( callRegionId );
        callRegion.subRegionCall( spentTime );

        if ( mCallTreeTable.get() != NULL )
        {
            int scl = 0; // not used
            mCallTreeTable->addCallCosts( callRegionId, regionId, scl, costs );
        }
    }
}

int TraceData::getCurrentRegionId( const char* regionName )
{
    int regionId = 0;  // default value avoids compiler warning due to exception

    if ( !mCallStack.empty() )
    {
        // check that regionName is the last region on the current call stack
        const int currentRegionId = mCallStack.currentRegionId();
        const RegionEntry& callRegion = getRegion( currentRegionId );

        if ( strcmp( callRegion.getRegionName(), regionName ) != 0 )
        {
            COMMON_THROWEXCEPTION( "mismatch in call stack, stop " << regionName <<
                                   ", but currently called region is " << callRegion.getRegionName() )
        }

        regionId = currentRegionId;
    }
    else
    {
        regionId = mRegionTable.getRegionId( regionName );
    }

    return regionId;
}

TraceData::TraceData( const char* prefix, ThreadId threadId, bool mThreadEnabled, bool callTreeFlag ) :
    mThreadId( threadId ),
    mRegionTable( mThreadEnabled ? common::thread::getThreadName( threadId )->c_str() : NULL )
{
    // calltree table only allocated if needed
    if ( callTreeFlag )
    {
        const char* threadName = mThreadEnabled ? common::thread::getThreadName( threadId )->c_str() : NULL;
        mCallTreeTable.reset( new CallTreeTable( prefix, threadName ) );
    }

    SCAI_LOG_DEBUG( logger, "TraceData for thread " << threadId )
}

TraceData::~TraceData()
{
    SCAI_LOG_DEBUG( logger, "~TraceData for thread " << mThreadId )
}

} /* end namespace tracing */

} /* end namespace scai */
