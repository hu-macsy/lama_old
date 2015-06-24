
#include "tracing/CallTree.hpp"
#include "tracing/CallStack.hpp"
#include "tracing/CallTreeTable.hpp"

#include "common/Exception.hpp"
#include "common/Walltime.hpp"

#include <fstream>
#include <stdint.h>
#include <unistd.h>
#include <vector>
#include <utility>

namespace tracing

{

using std::endl;

CallTreeTable CT;

LAMA_LOG_DEF_LOGGER( CallTree::logger, "CallTree" )

CallStack myCallStack;

void CallTree::enter( const int regionId, RegionEntry& region, const CounterArray& startValues )
{
    LAMA_LOG_DEBUG( logger, "enter " << regionId << " " << region.getRegionName() )

    if ( region.firstAccess() )
    {
        // write Info about the region in the file

        CT.writeRegion( regionId, region ) ;
    }

    if ( myCallStack.empty() )
    {
        // This is the first time we enter a region, do initialization

        CT.initCounters( startValues );
    }
    else
    {
        // Before we enter the called region add exclusive costs so far for caller region

        int caller_region = myCallStack.currentRegionId();

        int scl = 0;

        CT.addExclusiveCosts( caller_region, scl, startValues );
    }

    myCallStack.push( regionId, startValues );
}

void CallTree::leave( const int regionId, const RegionEntry& region, const CounterArray& counterVals )
{
    LAMA_LOG_DEBUG( logger, "leave " << regionId << " " << region.getRegionName() )

    // Add exlusive time for the current region 

    int scl = 0;

    CT.addExclusiveCosts( regionId, scl, counterVals );

    COMMON_ASSERT_EQUAL( myCallStack.currentRegionId(), regionId, "Call stack mismatch" )

    // Gets costs of the call and add this in CallTree table.

    CounterArray costs;

    myCallStack.getCosts( costs, counterVals );  // costs = counterVals - startVals

    myCallStack.pop();

    if ( !myCallStack.empty() )
    {
        CT.addCallCosts( myCallStack.currentRegionId(), regionId, scl, costs );
    }
}

}
