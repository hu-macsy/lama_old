/**
 * @file RegionTable.hpp
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
 * @brief Definition of class that contains all regions with their timings.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <common/config.hpp>
#include <common/Thread.hpp>

// others
#include <tracing/RegionEntry.hpp>
#include <tracing/CallStack.hpp>

// logging
#include <logging.hpp>

#include <cstring>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>

namespace tracing
{

/** Class that manages all regions in a table */

class COMMON_DLL_IMPORTEXPORT RegionTable
{

public:

    /** Constructor of a new region table.
     *
     *  The region table must contain the thread id as it might be
     *  written later by main thread.
     */

    RegionTable( const char* threadName );

    /** Destructor. */

    ~RegionTable();

    /** Get the id of a region, creates a new entry if region is not available yet.
     *
     *  @param[in] name   is the name of the region
     *  @param[in] file   is the name of the source file where region is defined
     *  @param[in] lno    is the line number of the region in the source file
     */

    int getRegionId( const char* name, const char* file, int lno );

    /** Get the id of an existing region, throws exception if not defined. */

    int getRegionId( const char* name );

    /** Return the elapsed time up to now.
     *  It will add also the time of a running region.
     */

    double elapsed( int regionId );

    /** Get full region record by its id. */

    RegionEntry& getRegion( int regionId );

    /** Get full region record by its id. */

    const RegionEntry& getRegion( int regionId ) const;

    /** This routine prints timing information. */

    void printTimer( std::ostream& );

    /** This routine prints timing information in std::cout */

    void printTimer();

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    struct    CmpString
    {
        bool operator()( const char* a, const char* b ) const
        {
            return std::strcmp( a, b ) < 0;
        }
    };

    // CallStack callStack;

    std::vector<RegionEntry> array; //!<  Entries for all timers

    /** Map of region strings to region ids that are the indexes to array.
     *
     *  Timer strings are given by pointer; pointers will be always valid as string
     *  remains as member variable in array.
     */

    std::map<const char*, int, CmpString> mapTimer;

    std::string mThreadName;

};

}
