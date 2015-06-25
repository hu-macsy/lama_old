/**
 * @file CallTreeTable.hpp
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
 * @brief Definition of class that keeps call tree informations and writes it to a file.
 * @author Thomas Brandes
 * @date 24.06.2015
 */
#pragma once

#include "tracing/RegionEntry.hpp"
#include "tracing/FileTable.hpp"
#include "tracing/Counters.hpp"

#include "logging/logging.hpp"

namespace tracing

{

/** Structure defines entry for call tree table.
 *
 *  An entry can be exclusive costs for a region (callee == -1) or
 *  costs for a certain call (region caller calls region callee).
 *
 *  The scl (source code location) can be set to distinguish between
 *  different calls or parts in the source code. The value 0 should be
 *  used if no distinction is made.
 */

struct CTTEntry
{
    int caller; // region id of the calling region
    int callee; // region id of the called region
    int scl;    // source code location for costs or for the call
    int calls;  // number of calls, as entries might be accumulated

    CounterArray costs;  // costs of a call or exclusive costs of a region

    void writeEntry( std::ostream& outfile );

    static void writeRegion( std::ostream& outfile, const int regionId, const int fileId, const RegionEntry& region );

    void addCallCosts( const CounterArray& callCosts )
    {
        costs += callCosts;
        calls++;  // addCallCosts stands for exactly one additional call
    }

    void set( int caller, int callee, int scl, const CounterArray& vals )
    {
        this->caller = caller;
        this->callee = callee;
        this->scl    = scl;
        costs = vals;
        calls = 1;
    }

    bool isSame( int other_caller, int other_callee, int other_scl )
    {
        return other_caller == caller && other_callee == callee && other_scl == scl;
    }
};

#define CALL_CACHE_SIZE 16

/** A CallTreeTable does not directly write all info records in the file but keeps
 *  them for a certain time in a cache.

 *  Note: each thread has its own table, so no synchronization is required here
 */

class CallTreeTable : private common::NonCopyable
{
public:

    /** Generate a new CallTree table and open the corresponding output file.
     *
     *  @param threadName is the name of the thread
     *
     *  threadName == NULL might be used if call tree data is only collected for
     *  the master/main thread.
     */

    CallTreeTable( const char* threadName );

    /** Destructor will also write final entries and close the output file. */

    ~CallTreeTable();

    /** Write an info record about a region in the file. */

    void writeRegion( const int regionId, const RegionEntry& region );

    /** Generate an info record about exclusive costs for a region */

    void addExclusiveCosts( const int regionId, const int scl, const CounterArray& currentCounterValues );

    /** Generate a cost call record about call costs for a region
     *
     *  @param[in] caller is the id of the calling region
     *  @param[in] callee is the id of the called region
     *  @param[in] scl source code line to distinguish between different calls in the source code
     *  @param[in] costs are counter values from entry to exit of the called region.
     */

    void addCallCosts( int caller, int callee, int scl, const CounterArray& costs );

    /** Initialize counter values for lastest counter values; are used to get exclusive costs. */

    void initCounters( const CounterArray& counterValues );

protected:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

private:

    CounterArray lastCounterValues;
    CounterArray totalCosts;

    CTTEntry mCallEntryCache[CALL_CACHE_SIZE];

    int callCachePos;
    int callCacheLast;

    int cacheHit;    // count hits for reused entry
    int cacheMiss;   // count misses for every new entry

    std::ofstream outfile;

    std::string mFileName;  // name of the output file

    int newPos();

    int find( int caller, int callee, int scl );

    void add( int caller, int callee, int scl, const CounterArray& costs );

    /** Clear the cache by writing out all entries. */

    void clear();

    /** Close output file and write final information. */

    void close();

    /** Open the output file for the calltree */

    void open( const char* threadName );

    FileTable mFileTable;
};

} // namespace
