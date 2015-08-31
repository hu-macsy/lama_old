/**
 * @file RegionTable.cpp
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
 * @brief Implementation of methods for class RegionTable.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/RegionTable.hpp>

// others
#include <scai/tracing/VTInterface.hpp>

#include <scai/common/Walltime.hpp>
#include <scai/common/Thread.hpp>
#include <scai/common/exception/Exception.hpp>

#include <cstdio>
#include <cstdlib>

using namespace std;

using scai::common::Thread;

namespace scai
{

namespace tracing
{

/* -------------------------------------------------------------------------- */

static Thread::Mutex printMutex; // needed to avoid mixing output of threads

/* -------------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( RegionTable::logger, "RegionTable" )

/* ---------------------------------------------------------------------- */

RegionTable::RegionTable( const char* threadName )
{
    // threadName can be NULL if tracing is only for main thread
    if ( threadName != NULL )
    {
        mThreadName = threadName;
    }

    SCAI_LOG_DEBUG( logger, "Constructor RegionTable" )
    // avoid too much reallocations at the beginning
    array.reserve( 16 );
}

/* ---------------------------------------------------------------------- */

RegionTable::~RegionTable()
{
    SCAI_LOG_DEBUG( logger, "~RegionTable" )

    if ( SCAI_LOG_INFO_ON( logger ) )
    {
        printTimer();
    }
}

/* ---------------------------------------------------------------------- */

double RegionTable::elapsed( int regionId )
{
    double elapsedTime = array[regionId].getInclusiveTime();
    /* Not yet

    // might be that region is still on the stack

    for( size_t i = 0; i < callStack.size(); i++ )
    {
        CallEntry& call = callStack[i];

        if( call.mRegion == regionId )
        {
            elapsedTime += common::Walltime::get() - call.mTimeStart;
        }
    }
    */
    return elapsedTime;
}

/* ---------------------------------------------------------------------- */

int RegionTable::getRegionId( const char* regionName, const char* file, int scl )
{
    std::map<const char*, int, CmpString>::iterator it = mapTimer.find( regionName );

    if ( it == mapTimer.end() )
    {
        int regionId = array.size();
        array.push_back( RegionEntry() );
        RegionEntry& entry = array[regionId];
        entry.mName = regionName;
        entry.mFile = file;
        entry.mLine = scl;
        VTInterface::define( entry );
        // do not use this: mapTimer[regionName] = regionId; // causes problems with composed strings
        mapTimer.insert( std::pair<const char*, int>( entry.mName.c_str(), regionId ));
        SCAI_LOG_DEBUG( logger,
                        "Added region " << regionId << "( " << array[regionId].mName << ") for region " << regionName )
        return static_cast<int>( regionId );
    }
    else
    {
        // alread defined
        return it->second;
    }
}

/* ---------------------------------------------------------------------- */

int RegionTable::getRegionId( const char* regionName )
{
    std::map<const char*, int, CmpString>::iterator it = mapTimer.find( regionName );

    if ( it == mapTimer.end() )
    {
        COMMON_THROWEXCEPTION( "Region " << regionName << " never defined" )
    }
    else
    {
        // already defined
        return it->second;
    }
}

/* ---------------------------------------------------------------------- */

RegionEntry& RegionTable::getRegion( int regionId )
{
    return array[regionId];
}

const RegionEntry& RegionTable::getRegion( int regionId ) const
{
    return array[regionId];
}

/* ---------------------------------------------------------------------- */

void RegionTable::printTimer()
{
    Thread::ScopedLock lock( printMutex );
    printTimer( cout );
}

/* ---------------------------------------------------------------------- */

void RegionTable::printTimer( ostream& outfile )
{
    SCAI_LOG_INFO( logger, "Summary of all timers for thread " << mThreadName  );

    if ( mThreadName.size() > 0 )
    {
        outfile << "========================================" << endl;
        outfile << "Timing info of regions for Thread " << mThreadName << endl;
        outfile << "========================================" << endl;
    }
    else
    {
        outfile << "=========================================" << endl;
        outfile << "Timing info of regions (only main thread)" << endl;
        outfile << "=========================================" << endl;
    }

    // use map iterator for alphabetical output
    std::map<const char*, int, CmpString>::iterator it;

    for ( it = mapTimer.begin(); it != mapTimer.end(); it++ )
    {
        int regionId = it->second;
        const RegionEntry& region = array[regionId];
        region.printTime( outfile );
    }
}

} /* end namespace tracing */

} /* end namespace scai */
