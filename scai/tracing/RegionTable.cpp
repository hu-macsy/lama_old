/**
 * @file RegionTable.cpp
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
 * @brief Implementation of methods for class RegionTable.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

// hpp
#include <scai/tracing/RegionTable.hpp>

// local library
#include <scai/tracing/VTInterface.hpp>

// internal scai libraries
#include <scai/common/Walltime.hpp>
#include <scai/common/macros/throw.hpp>

// std
#include <cstdio>
#include <cstdlib>

using namespace std;

namespace scai
{

namespace tracing
{

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
        printTimer( cout );
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
    MapRegion::iterator it = mapTimer.find( regionName );

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
        mapTimer.insert( std::pair<const char*, int>( entry.mName.c_str(), regionId ) );
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
    MapRegion::iterator it = mapTimer.find( regionName );

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
    MapRegion::iterator it;

    for ( it = mapTimer.begin(); it != mapTimer.end(); it++ )
    {
        int regionId = it->second;
        const RegionEntry& region = array[regionId];
        region.printTime( outfile );
    }
}

} /* end namespace tracing */

} /* end namespace scai */
