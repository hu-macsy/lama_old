/**
 * @file CallTreeTable.cpp
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
 * @brief Implementation of class that keeps call tree informations and writes it to a file.
 * @author Thomas Brandes
 * @date 24.06.2015
 */

#include "common/Exception.hpp"

#include "tracing/CallTreeTable.hpp"

namespace tracing
{

using namespace std;

LAMA_LOG_DEF_LOGGER( CallTreeTable::logger, "CallTreeTable" )

void CTTEntry::writeEntry( ostream& outfile )
{
    if ( callee == - 1 )
    {
        // cost line
        outfile << "# begin exclusive cost line" << endl;
        outfile << "fl 0" << endl;
        outfile << "fn " << caller << endl;
        outfile << scl << " ";
        costs.write( outfile, " " );
        outfile << endl;
        outfile << "# end exclusive cost line" << endl;
    }
    else
    {
        outfile << "# begin call cost line" << endl;
        outfile << "fl 0" << endl;
        outfile << "fn " << caller << endl;
        outfile << "cfl 0" << endl;
        outfile << "cfn " << callee << endl;
        outfile << "calls " << calls << " 0" << endl;
        outfile << scl << " ";
        costs.write( outfile, " " );
        outfile << endl;
        outfile << "# end call cost line" << endl;
    }
}

void CTTEntry::writeRegion( ostream& outfile, const int regionId, const RegionEntry& region )
{
    CounterArray zeroCosts;   // all values are zero
    outfile << "# begin region info line" << endl;
    outfile << "fl " << region.getFileToken() << " " << region.getFileName() << endl;
    outfile << "fn " << regionId;
    outfile << " 0";  //  src_regionId, does not matter
    outfile << " " <<  region.getRegionName();
    outfile << " 2";  // function
    outfile << " " << region.getLine();
    outfile << " " << region.getLine();
    outfile << " ?";   // module name, does not matter
    outfile << endl;
    outfile << "0 ";
    zeroCosts.write( outfile, " " );
    outfile << endl;
    outfile << "# end region info line" << endl;
}

int CallTreeTable::newPos()
{
    if ( callCachePos < CALL_CACHE_SIZE )
    {
        return callCachePos++;
    }

    call_cache[callCacheLast].writeEntry( outfile );
    int pos = callCacheLast++;

    if ( callCacheLast == CALL_CACHE_SIZE )
    {
        callCacheLast = 0;
    }

    return pos;
}

/** Write out outstanding entries. */

void CallTreeTable::clear()
{
    for ( int i = 0; i < callCachePos; ++i )
    {
        call_cache[i].writeEntry( outfile );
    }

    callCachePos  = 0;
    callCacheLast = 0;
}

void CallTreeTable::close()
{
    // file might not have been opened, can be closed several times

    if ( outfile.is_open() )
    {
        outfile << "# closed by close" << endl;
        outfile << "totals ";
        totalCosts.write( outfile, " " );
        outfile << endl;
        outfile << "# cache has " << cacheHit << " hits and "
                << cacheMiss << " misses" << endl;
        outfile.close();
    }
}

void CallTreeTable::open()
{
    if ( outfile.is_open() )
    {
        return;
    }

    LAMA_LOG_DEBUG( logger, "open calltree file" );
    outfile.open( "calltree.ct", ios::out );

    if ( outfile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open calltree.ct" )
    }

    outfile << "pid " << getpid() << endl;
    outfile << "cmd program_name" << endl;
    outfile << "part 1" << endl;
    outfile << endl;
    outfile << "# rate " << common::Walltime::timerate() << endl;
    double rate = 1.0 / common::Walltime::timerate();
    outfile << "event WALL_TICKS wallticks" << endl;
    outfile << "events WALL_TICKS" << endl;
    outfile << "define WALL_TIME " << rate << " WALL_TICKS" << endl;
}

CallTreeTable::CallTreeTable()
{
    callCachePos = 0;
    callCacheLast = 0;
    cacheHit = 0;
    cacheMiss = 0;
    open();
}

CallTreeTable::~CallTreeTable()
{
    clear();   // write all entries not written yet
    close();   // close the outfile
}

void CallTreeTable::writeRegion( const int regionId, const RegionEntry& region )
{
    CTTEntry::writeRegion( outfile, regionId, region );
}

int CallTreeTable::find( int caller, int callee, int scl )
{
    for ( int i = 0; i < callCachePos; ++i )
    {
        if ( call_cache[i].isSame( caller, callee, scl ) )
        {
            cacheHit++;
            return i;
        }
    }

    cacheMiss++;
    return -1;
}

void CallTreeTable::add( int caller, int callee, int scl, const CounterArray& costs )
{
    int pos = find( caller, callee, scl );

    if ( pos >= 0 )
    {
        call_cache[pos].addCallCosts( costs );
    }
    else
    {
        pos = newPos();
        call_cache[pos].set( caller, callee, scl, costs );
    }
}

void CallTreeTable::addExclusiveCosts( const int regionId, const int scl, const CounterArray& currentCounterValues )
{
    CounterArray costs = currentCounterValues - lastCounterValues;

    totalCosts += costs;     // sum of all exclusive costs will be the total costs

    lastCounterValues = currentCounterValues;  // uses current counters as start for next time

    add( regionId, -1, scl, costs );
}

void CallTreeTable::addCallCosts( int caller, int callee, int scl, const CounterArray& costs )
{
    add( caller, callee, scl, costs );
}

void CallTreeTable::initCounters( const CounterArray& counterValues )
{
    lastCounterValues = counterValues;
}

} // namespace
