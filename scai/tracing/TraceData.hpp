/**
 * @file TraceData.hpp
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
 * @brief Definition of class that contains all data used for tracing.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

#include <scai/tracing/RegionTable.hpp>
#include <scai/tracing/CallStack.hpp>
#include <scai/tracing/CallTreeTable.hpp>

namespace scai
{

namespace tracing
{

/** Structure that contains all relevant data structures used for tracing. */

class TraceData
{
public:

    typedef common::Thread::Id ThreadId;

    /** Constructor of new thread record for tracing data.
     *
     *  @param[in] prefix string for the first part of the calltree filename
     *  @param[in] threadId  id of the thread to which data belongs
     *  @param[in] threadEnabled if true tracing is done for all threads
     */

    TraceData( const char* prefix, ThreadId threadId, bool threadEnabled );

    /** Destructor. */

    ~TraceData();

    /** Query of the thread id to which thread data belongs. */

    ThreadId getId() const
    {
        return mThreadId;
    }

    void leave( const int regionId, RegionEntry& region, const bool callTreeFlag );

    void enter( const int regionId, RegionEntry& region, const bool callTreeFlag );

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

    CallTreeTable mCallTreeTable;   // generates call tree trace file

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

};

} /* end namespace tracing */

} /* end namespace scai */
