
#pragma once

#include "tracing/RegionEntry.hpp"
#include "tracing/Counters.hpp"

#include "logging/logging.hpp"

#include <vector>

namespace tracing

{

class CallTree 
{
public:

    static void enter( const int region_id, RegionEntry& region, const CounterArray& startValues );

    static void leave( const int region_id, const RegionEntry& region, const CounterArray& stopValues );

    static void finish();

private:

    static std::vector<int> theCallStack;

    LAMA_LOG_DECL_STATIC_LOGGER( logger )
};

}
