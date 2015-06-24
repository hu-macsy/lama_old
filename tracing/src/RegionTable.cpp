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
#include <tracing/RegionTable.hpp>

// others
#include <tracing/VTInterface.hpp>

#include <common/Walltime.hpp>
#include <common/Thread.hpp>
#include <common/Exception.hpp>

#include <cstdio>
#include <cstdlib>

using namespace std;

using common::Thread;

namespace tracing
{

/* -------------------------------------------------------------------------- */

static Thread::Mutex printMutex; // needed to avoid mixing output of threads

/* -------------------------------------------------------------------------- */

LAMA_LOG_DEF_LOGGER( RegionTable::logger, "RegionTable" )

/* ---------------------------------------------------------------------- */

RegionTable::RegionTable( ThreadId threadId )
    : mThreadId( threadId )
{
    LAMA_LOG_DEBUG( logger, "Thread " << threadId << ": creates RegionTable" )
}

/* ---------------------------------------------------------------------- */

RegionTable::~RegionTable()
{
    if( LAMA_LOG_INFO_ON( logger ) )
    {
        printTimer();
    }
}

/* ---------------------------------------------------------------------- */

void RegionTable::init()
{
    // clear the callstack

    callStack.clear();
}

/* ---------------------------------------------------------------------- */

void RegionTable::start( int regionId, const CounterArray& startValues )
{
    LAMA_LOG_DEBUG( logger, "start region " << regionId << ", time = " << startValues.getWalltime() )

    callStack.push( regionId, startValues );
}

/* ---------------------------------------------------------------------- */

void RegionTable::stop( int regionId, const CounterArray& stopValues )
{
    LAMA_LOG_DEBUG( logger, "stop region " << regionId << ", time = " << stopValues.getWalltime() )

    if( callStack.empty() )
    {
        LAMA_LOG_ERROR( logger, "stop region on empty call region stack" )
        return;
    }

    const int currentRegionId = callStack.currentRegionId();

    if( ( currentRegionId >= 0 ) )
    {
        COMMON_ASSERT_EQUAL( currentRegionId, regionId, 
                             "mismatch call stack, current region = " << array[currentRegionId].mName
                             << ", stop for " << array[regionId].mName )
    }

    double spentTime = stopValues.getWalltime( callStack.currentCounters() );

    array[regionId].addCall( spentTime );

    callStack.pop();

    // correct exclusive time of previous entry in call stack

    if( !callStack.empty() )
    {
        int currentRegionId = callStack.currentRegionId();
        array[currentRegionId].subRegionCall( spentTime );
    }
}

/* ---------------------------------------------------------------------- */

int RegionTable::getCurrentRegionId( const char* regionName )
{
    if( !callStack.empty() )
    {
        // check that regionName is the last region on the current call stack

        const int currentRegionId = callStack.currentRegionId();

        const RegionEntry& callRegion = getRegion( currentRegionId );

        if( strcmp( callRegion.getRegionName(), regionName ) != 0 )
        {
            fprintf( stderr, "mismatch in call stack, stop %s but in %s\n", regionName, callRegion.getRegionName() );
            return 0;
        }

        return currentRegionId;
    }
    else
    {
        // no callstack available, so we search it

        std::map<const char*,int,CmpString>::iterator it = mapTimer.find( regionName );

        if( it == mapTimer.end() )
        {
            LAMA_LOG_FATAL( logger, "Region " << regionName << " never defined" )
            exit( -1 );
        }
        else
        {
            // alread defined

            return it->second;
        }
    }

    return 0; // avoids warning
}

/* ---------------------------------------------------------------------- */

void RegionTable::stop( const char* regionName, const CounterArray& counterValues )
{
    // check that regionName is the last region on the current call stack

    const int currentRegionId = callStack.currentRegionId();

    const RegionEntry callRegion = getRegion( currentRegionId );

    if( strcmp( callRegion.getRegionName(), regionName ) != 0 )
    {
        fprintf( stderr, "mismatch in call stack, stop %s but in %s\n", regionName, callRegion.getRegionName() );
        return;
    }

    stop( currentRegionId, counterValues );
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

int RegionTable::getRegion( const char* id, const char* file, int lno )
{
    std::map<const char*,int,CmpString>::iterator it = mapTimer.find( id );

    if( it == mapTimer.end() )
    {
        const size_t regionId = array.size();
        array.resize( regionId + 1 );
        RegionEntry& entry = array[regionId];
        entry.mName = id;
        entry.mFile = file;
        entry.mLine = lno;

        VTInterface::define( entry );

        // do not use this: mapTimer[id] = regionId; // causes problems with composed strings
        mapTimer[entry.mName.c_str()] = static_cast<int>( regionId );
        LAMA_LOG_DEBUG( logger,
                        "Thread " << mThreadId << " added region " << regionId << "( " << array[regionId].mName << ") for region " << id )
        return static_cast<int>( regionId );
    }
    else
    {
        // alread defined

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
    LAMA_LOG_INFO( logger, "Summary of all timers for thread " << Thread::getThreadId( mThreadId ) );

    outfile << "========================================" << endl;
    outfile << "Timing info of regions for Thread " << Thread::getThreadId( mThreadId ) << endl;
    outfile << "========================================" << endl;

    // use map iterator for alphabetical output

    std::map<const char*,int,CmpString>::iterator it;

    for( it = mapTimer.begin(); it != mapTimer.end(); it++ )
    {
        int regionId = it->second;
        const std::string& name = it->first;

        const RegionEntry& region = array[regionId];

        region.printTime( outfile );
    }
}

} // namespace tracing
