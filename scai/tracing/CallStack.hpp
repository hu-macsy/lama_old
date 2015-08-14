/**
 * @file CallStack.hpp
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
 * @brief Definition of class for CallStack
 * @author Thomas Brandes
 * @date 23.06.2015
 */

#pragma once

#include <scai/tracing/Counters.hpp>
#include <vector>

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

}  // namespace
