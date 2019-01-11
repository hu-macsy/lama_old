/**
 * @file RegionTable.hpp
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
 * @brief Definition of class that contains all regions with their timings.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai library
#include <scai/tracing/RegionEntry.hpp>
#include <scai/tracing/CallStack.hpp>

#include <scai/logging.hpp>

// std
#include <cstring>
#include <cstdio>
#include <vector>
#include <map>
#include <fstream>

namespace scai
{

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

private:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

    // CallStack callStack;

    std::vector<RegionEntry> array; //!<  Entries for all timers

    typedef std::map<std::string, int> MapRegion;

    /** Map of region strings to region ids that are the indexes to array.
     *
     *  Timer strings are stored as std::string and not as const char*
     */

    MapRegion mapTimer;

    std::string mThreadName;
};

} /* end namespace tracing */

} /* end namespace scai */
