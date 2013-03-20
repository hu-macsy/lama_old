/**
 * @file RegionTable.hpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Definition of class that contains all regions with their timings.
 * @author Thomas Brandes
 * @date 21.11.2011
 * $Id$
 */

#ifndef LAMA_REGION_TABLE_HPP_
#define LAMA_REGION_TABLE_HPP_

// for dll_import
#include <lama/config.hpp>

#include <lama/tracing/RegionEntry.hpp>

// others
#include <lama/task/Thread.hpp>

// logging
#include <logging/logging.hpp>

#include <cstring>
#include <cstdio>
#include <vector>
#include <map>

#define NEW_VT

namespace tracing
{

/** Class to collect timing information about different regions in
 an instrumented program.

 For each region timers are defined. Timing for a routine must
 be started and stopped according a subroutine call structure.

 \code
 start(region1)
 start(region2)
 stop(region2)
 start(region3)
 elapsed(region1)
 stop(region3)
 stop(region1)
 \endcode

 At the end of a run, inclusive and exclusive time is available for
 all regions.
 */

class LAMA_DLL_IMPORTEXPORT RegionTable
{

public:

    /** Constructor of a new region table.
     *
     *  @param[in] threadId  id of the thread to which region table belongs
     */

    RegionTable( lama::Thread::Id threadId );

    /** Destructor. */

    ~RegionTable();

    /** Query of the thread id to which thread table belongs. */

    lama::Thread::Id getId() const
    {
        return mThreadId;
    }

    /** Get the id of a region, creates a new entry if region is not available yet.
     *
     *  @param[in] id     is the name of the region
     *  @param[in] file   is the name of the source file where region is defined
     *  @param[in] lno    is the line number of the region in the source file
     */

    int getRegion( const char* id, const char* file, int lno );

    /** Get the id of the current region
     *
     *  @param[in] regionName   must match the name of the current region
     */
    int getCurrentRegionId( const char* regionName );

    /** Initialization of the timers. Resets all timers. */

    void init();

    /** Enter a region with the timestamp. */

    void start( int regionId, double wallTime );

    /** Leave the region with the timestamp.
     *
     *  @param[in] regionId is the id of region to leave, must match the last on call stack
     *  @param[in] wallTime is time stamp of call
     */
    void stop( int regionId, double wallTime );

    /** Leave the region with the timestamp. */

    void stop( const char* regionName, double wallTime );

    /** Return the elapsed time up to now.
     It will add also the time of a running region.
     */

    double elapsed( int regionId );

    /** Get full region record by its id. */

    RegionEntry& getRegion( int regionId );

    /** Get full region record by its id. */

    const RegionEntry& getRegion( int regionId ) const;

    /** This routine prints timing information. */

    void printTimer();

    void printTimer( FILE* f );

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger );

    struct CmpString
    {
        bool operator()( const char* a, const char* b ) const
        {
            return std::strcmp( a, b) < 0;
        }
    };

    /** Structure for entry on the call stack (id, starttime) */

    struct CallEntry
    {
        CallEntry( int region, double timeStart )
        {
            mRegion = region;
            mTimeStart = timeStart;
        }

        int mRegion; //!< reference id for the region
        double mTimeStart;//!< absolute time when call started
    };

    std::vector<CallEntry> callStack;

    std::vector<RegionEntry> array; //!<  Entries for all timers

    lama::Thread::Id mThreadId;

    /** Map of region strings to region ids that are the indexes to array.
     *
     *  Timer strings are given by pointer; pointers will be always valid as string
     *  remains as member variable in array.
     */

    std::map<const char*, int, CmpString> mapTimer;
};

}

#endif //LAMA_REGION_TABLE_HPP_
