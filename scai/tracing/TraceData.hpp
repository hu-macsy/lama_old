/**
 * @file TraceData.hpp
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
 * @brief Definition of class that contains all data used for tracing.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// local library
#include <scai/tracing/RegionTable.hpp>
#include <scai/tracing/CallStack.hpp>
#include <scai/tracing/CallTreeTable.hpp>

#include <scai/common/thread.hpp>

#include <memory>

namespace scai
{

namespace tracing
{

/** Structure that contains all relevant data structures used for tracing. */

class TraceData
{
public:

    typedef common::thread::Id ThreadId;

    /** Constructor of new thread record for tracing data.
     *
     *  @param[in] prefix string for the first part of the calltree filename
     *  @param[in] threadId  id of the thread to which data belongs
     *  @param[in] threadEnabled if true tracing is done for all threads
     *  @param[in] callTreeFlag if true call tree trace file will be generated
     */

    TraceData( const char* prefix, ThreadId threadId, bool threadEnabled, bool callTreeFlag );

    /** Destructor. */

    ~TraceData();

    /** Query of the thread id to which thread data belongs. */

    ThreadId getId() const
    {
        return mThreadId;
    }

    /** Enter a region, will set corresponding tracing values.
     *
     *  @param [in] regionId unique identification id of the region
     *  @param [in,out] region is the record with trace data of the region
     */
    void enter( const int regionId, RegionEntry& region );

    /** Leave a region, will set corresponding tracing values.
     *
     *  @param [in] regionId unique identification id of the region
     *  @param [in,out] region is the record with trace data of the region
     */
    void leave( const int regionId, RegionEntry& region );

    /** Get the id of a region, creates a new entry if region is not available yet.
     *
     *  @param[in] regionName is the name of the region
     *  @param[in] file       is the name of the source file where region is defined
     *  @param[in] scl        is the line number of the region in the source file
     */
    int getRegionId( const char* regionName, const char* file, int scl )
    {
        return mRegionTable.getRegionId( regionName, file, scl );
    }

    /** Get the id of the current region
     *
     *  @param[in] regionName   must match the name of the current region
     */
    int getCurrentRegionId( const char* regionName );

    RegionEntry& getRegion( int regionId )
    {
        return mRegionTable.getRegion( regionId );
    }

    void printTimer( std::ostream& outfile )
    {
        mRegionTable.printTimer( outfile );
    }

private:

    ThreadId      mThreadId;

    CallStack     mCallStack;

    RegionTable   mRegionTable;     // counts timing of regions

    std::unique_ptr<CallTreeTable> mCallTreeTable;   // generates call tree trace file

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

} /* end namespace tracing */

} /* end namespace scai */
