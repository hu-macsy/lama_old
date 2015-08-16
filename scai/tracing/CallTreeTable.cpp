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

#include <scai/tracing/CallTreeTable.hpp>

#include <scai/common/Exception.hpp>

// getpid required
#include <unistd.h>

namespace scai
{

namespace tracing
{

using namespace std;

SCAI_LOG_DEF_LOGGER( CallTreeTable::logger, "CallTreeTable" )

void CTTEntry::writeEntry( ostream& outfile )
{
    if ( callee == - 1 )
    {
        // cost entry
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
        // call entry
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

void CTTEntry::writeRegion( ostream& outfile, const int regionId, const int fileId, const RegionEntry& region )
{
    std::string regionName = region.getRegionName();
    std::string groupName  = "?";
    size_t pindex = regionName.find_first_of( "." );

    if ( pindex != std::string::npos )
    {
        // split the region name and define it
        groupName  = regionName.substr( 0, pindex );
        regionName = regionName.substr( pindex + 1 );
    }

    CounterArray zeroCosts;   // all values are zero
    outfile << "# begin region info line" << endl;
    outfile << "fl " << fileId << " " << region.getFileName() << endl;
    outfile << "fn " << regionId;
    outfile << " 0";  //  src_regionId, does not matter
    outfile << " " <<  regionName;
    outfile << " 2";  // function
    outfile << " " << region.getLine();
    outfile << " " << region.getLine();
    outfile << " " << groupName;
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

    mCallEntryCache[callCacheLast].writeEntry( outfile );
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
        mCallEntryCache[i].writeEntry( outfile );
    }

    callCachePos  = 0;
    callCacheLast = 0;
}

void CallTreeTable::close()
{
    // file might not have been opened, can be closed several times
    if ( outfile.is_open() )
    {
        SCAI_LOG_DEBUG( logger, "close calltree file " << mFileName );
        outfile << "# closed by close" << endl;
        outfile << "totals ";
        totalCosts.write( outfile, " " );
        outfile << endl;
        outfile << "# cache has " << cacheHit << " hits and "
                << cacheMiss << " misses" << endl;
        outfile.close();
        mFileName.clear();
    }
}

void CallTreeTable::open( const char* threadSuffix )
{
    if ( outfile.is_open() )
    {
        return;
    }

    mFileName = "calltree.ct";
    SCAI_LOG_DEBUG( logger, "open calltree file" );

    if ( threadSuffix != NULL )
    {
        mFileName += ".";
        mFileName += threadSuffix;
    }

    outfile.open( mFileName.c_str(), ios::out );

    if ( outfile.fail() )
    {
        COMMON_THROWEXCEPTION( "Could not open " << mFileName )
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

CallTreeTable::CallTreeTable( const char* threadSuffix )
{
    callCachePos = 0;
    callCacheLast = 0;
    cacheHit = 0;
    cacheMiss = 0;
    open( threadSuffix );
}

CallTreeTable::~CallTreeTable()
{
    clear();   // write all entries not written yet
    close();   // close the outfile
}

void CallTreeTable::writeRegion( const int regionId, const RegionEntry& region )
{
    int fileId = mFileTable.getFileId( region.getFileName() );
    CTTEntry::writeRegion( outfile, regionId, fileId, region );
}

int CallTreeTable::find( int caller, int callee, int scl )
{
    for ( int i = 0; i < callCachePos; ++i )
    {
        if ( mCallEntryCache[i].isSame( caller, callee, scl ) )
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
        mCallEntryCache[pos].addCallCosts( costs );
    }
    else
    {
        pos = newPos();
        mCallEntryCache[pos].set( caller, callee, scl, costs );
    }
}

void CallTreeTable::addExclusiveCosts( const int regionId, const int scl, const CounterArray& currentCounterValues )
{
    CounterArray costs = currentCounterValues - lastCounterValues;
    SCAI_LOG_DEBUG( logger, "region " << regionId << ", add exclusive costs: " << costs )
    totalCosts += costs;     // sum of all exclusive costs will be the total costs
    lastCounterValues = currentCounterValues;  // uses current counters as start for next time
    add( regionId, -1, scl, costs );
}

void CallTreeTable::addCallCosts( int caller, int callee, int scl, const CounterArray& costs )
{
    SCAI_LOG_DEBUG( logger, "call " << caller << " -> " << callee << ", add call costs: " << costs )
    add( caller, callee, scl, costs );
}

void CallTreeTable::initCounters( const CounterArray& counterValues )
{
    lastCounterValues = counterValues;
}

} /* end namespace tracing */

} /* end namespace scai */
