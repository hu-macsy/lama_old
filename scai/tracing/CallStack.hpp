/**
 * @file CallStack.hpp
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
 * @brief Definition of class for CallStack
 * @author Thomas Brandes
 * @date 23.06.2015
 */

#pragma once

// local library
#include <scai/tracing/Counters.hpp>

// std
#include <vector>

namespace scai
{

namespace tracing
{

/** Structure for an entry in the call stack. */

struct CallStackEntry
{
    int          mRegionId;       //! name of the region
    CounterArray mCounterValues;  //! values of performance counters at entry

    CallStackEntry( int regionId, const CounterArray& counterValues ) :
        mRegionId( regionId ),
        mCounterValues( counterValues )
    {
    }

    // Note: default copy constructor, assignment operator might be used.
};

/** The class CallStack contains the current nesting of calls. */

class CallStack
{
public:

    void clear()
    {
        stack.clear();
    }

    void push( int regionId, const CounterArray& entryCounterValues )
    {
        stack.push_back( CallStackEntry( regionId, entryCounterValues ) );
    }

    void pop()
    {
        stack.pop_back();
    }

    /** Get the costs for the total call of the current region
     *
     *  @param[out] costs are the difference between counterValues at entry and counterVals
     *  @param[in]  counterVals are counter values at exit of the current region
     */

    void getCosts( CounterArray& costs, const CounterArray& counterVals ) const
    {
        const CounterArray& enterVals = stack.back().mCounterValues;
        costs = counterVals - enterVals;
    }

    /** Get the current region id, taken from last entry of call stack */

    int currentRegionId() const
    {
        return stack.back().mRegionId;
    }

    const CounterArray& currentCounters() const
    {
        return stack.back().mCounterValues;
    }

    /** Query if call stack is empty. */

    bool empty() const
    {
        return stack.empty();
    }

private:

    std::vector<CallStackEntry> stack;
};

} /* end namespace tracing */

} /* end namespace scai */
