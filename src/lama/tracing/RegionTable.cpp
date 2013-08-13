/**
 * @file RegionTable.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @date 21.11.2011
 * @since 1.0.0
 */

// hpp
#include <lama/tracing/RegionTable.hpp>

// others
#include <lama/tracing/VTInterface.hpp>
#include <lama/Walltime.hpp>
#include <lama/exception/LAMAAssert.hpp>

#include <cstdio>

namespace tracing
{

/* -------------------------------------------------------------------------- */

static boost::mutex printMutex; // needed to avoid mixing output of threads

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
    if LAMA_LOG_INFO_ON( logger )
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

void RegionTable::start( int regionId, double wallTime )
{
    LAMA_LOG_DEBUG( logger, "start region " << regionId << ", time = " << wallTime )

    callStack.push_back( CallEntry( regionId, wallTime ) );
}

/* ---------------------------------------------------------------------- */

void RegionTable::stop( int regionId, double wallTime )
{
    LAMA_LOG_DEBUG( logger, "stop region " << regionId << ", time = " << wallTime )

    if ( callStack.size() == 0 )
    {
        LAMA_LOG_ERROR( logger, "stop region on empty call region stack" )
        return;
    }

    CallEntry call = callStack.back();

    if ( ( call.mRegion >= 0 ) && ( call.mRegion != regionId ) )
    {
        fprintf( stderr, "mismatch in call stack, stop %s but in %s\n", array[regionId].mName.c_str(),
                 array[call.mRegion].mName.c_str() );
        return;
    }

    callStack.pop_back();

    double spentTime = wallTime - call.mTimeStart;

    array[regionId].addCall( spentTime );

    // correct exclusive time of previous entry in call stack

    if ( callStack.size() )
    {
        CallEntry before = callStack.back();
        array[before.mRegion].subRegionCall( spentTime );
    }
}

/* ---------------------------------------------------------------------- */

int RegionTable::getCurrentRegionId( const char* regionName )
{
    if ( callStack.size() )
    {
        // check that regionName is the last region on the current call stack

        const CallEntry& call = callStack.back();

        const RegionEntry& callRegion = getRegion( call.mRegion );

        if ( strcmp( callRegion.getRegionName(), regionName ) != 0 )
        {
            fprintf( stderr, "mismatch in call stack, stop %s but in %s\n", regionName, callRegion.getRegionName() );
            return 0;
        }

        return call.mRegion;
    }
    else
    {
         // no callstack available, so we search it

        std::map<const char*,int,CmpString>::iterator it = mapTimer.find( regionName );

        if ( it == mapTimer.end() )
        {
            LAMA_THROWEXCEPTION( "Region " << regionName << " never defined" )
        }
        else
        {
            // alread defined

            return it->second;
        }
    }

    return 0;  // avoids warning
}

/* ---------------------------------------------------------------------- */

void RegionTable::stop( const char* regionName, double wallTime )
{
    // check that regionName is the last region on the current call stack

    const CallEntry& call = callStack.back();

    const RegionEntry callRegion = getRegion( call.mRegion );

    if ( strcmp( callRegion.getRegionName(), regionName ) != 0 )
    {
        fprintf( stderr, "mismatch in call stack, stop %s but in %s\n", regionName, callRegion.getRegionName() );
        return;
    }

    stop( call.mRegion, wallTime );
}

/* ---------------------------------------------------------------------- */

double RegionTable::elapsed( int regionId )
{
    double elapsedTime = array[regionId].getInclusiveTime();

    // might be that region is still on the stack

    for ( size_t i = 0; i < callStack.size(); i++ )
    {
        CallEntry& call = callStack[i];

        if ( call.mRegion == regionId )
        {
            elapsedTime += lama::Walltime::get() - call.mTimeStart;
        }
    }

    return elapsedTime;
}

/* ---------------------------------------------------------------------- */

int RegionTable::getRegion( const char* id, const char* file, int lno )
{
    std::map<const char*,int,CmpString>::iterator it = mapTimer.find( id );

    if ( it == mapTimer.end() )
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
    boost::mutex::scoped_lock scoped_lock( printMutex );

    printTimer( stdout );
}

/* ---------------------------------------------------------------------- */

void RegionTable::printTimer( FILE* f )
{
    std::cout << "Summary of all timers for thread " << mThreadId << std::endl;

    std::ostringstream threadInfo; 

    threadInfo << mThreadId;  
    
    fprintf( f, "====================================\n" );
    fprintf( f, "Timing info of regions for Thread %s\n", threadInfo.str().c_str() );
    fprintf( f, "====================================\n" );

    // use map iterator for alphabetical output

    std::map<const char*,int,CmpString>::iterator it;

    for ( it = mapTimer.begin(); it != mapTimer.end(); it++ )
    {
        int regionId = it->second;
        const std::string& name = it->first;

        const RegionEntry& region = array[regionId];

        fprintf( f, "Time %s (in ms) : #calls = %d, inclusive = %5.2f, exclusive = %5.2f\n", name.c_str(),
                 region.getCalls(), region.getInclusiveTime() * 1000.0, region.getExclusiveTime() * 1000.0 );
    }
}

} // namespace tracing
