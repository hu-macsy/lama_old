
#pragma once

#include "tracing/RegionEntry.hpp"

#include <vector>

namespace tracing

{

class CallTree 
{
public:

    static void enter( const int region_id, RegionEntry& region );

    static void leave( const int region_id, const RegionEntry& region );

    static void finish();

private:

    static std::vector<int> theCallStack;
};

}
