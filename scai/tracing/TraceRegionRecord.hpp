/**
 * @file TraceRegionRecord.hpp
 *
 * @license
 * Copyright (c) 2009-2017
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Definition of class for tracing that measures inclusive and exclusive
 *        walltimes for regions defined by directive SCAI_REGION.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>

// internal scai libraries
#include <scai/logging.hpp>

#include <memory>

namespace scai
{

/** Namespace for all data structures used in tracing library. */

namespace tracing
{

/** This class is a helper class for tracing a scope in C++.
 *
 *  The constructor of an object creates a start entry and the destructor the stop entry
 *  for trace files. The use of this class is much safer than using calls of routines start and
 *  stop explicitly as it works also very well in case of exception.
 */

class COMMON_DLL_IMPORTEXPORT TraceRegionRecord
{

public:

    /** Constructor of a tracer object, produces 'start' entry in trace file.
     *
     *  @param[in] regionName  name of the region, must be unique in application
     *  @param[in] fileName  name of the file in which region is coded
     *  @param[out] lno  line number in file where the region starts
     */

    TraceRegionRecord( const char* regionName, const char* fileName, int lno );

    /** Same as before but additional integer value that is used as suffix for region name.
     */

    TraceRegionRecord( const char* regionName, int n, const char* fileName, int lno );

    /** Trace record for the current region, name must match regionName */

    TraceRegionRecord( const char* regionName );

    /** Destructor of a tracer object, produces 'stop' entry in trace file. */

    ~TraceRegionRecord();

    /**
     *  Query for the inclusive time spent for the last call of a region.
     *
     *  @param[in] regionName name of the region that is queried
     *  @return inclusive time spent for the last call of region
     *
     *  This method must not be called within the region itself, i.e. region
     *  must not be on the current call stack.
     */
    static double spentLast( const char* regionName );

    /** Generate only a 'start region' entry in trace file. */

    static void start( const char* regionName, const char* file, int lno );

    /** Generate only a 'end region' entry in trace file. */

    static void stop( const char* regionName );

    void enter();

    void leave();

protected:

    SCAI_LOG_DECL_STATIC_LOGGER( logger )

private:

    /** Common routine for all constructors to check for settings of trace. */

    void initSettings();

    /** Each region timing keeps a shared pointer to the configuration.
     *  By this way it is guaranteed that timer information is only printed
     *  when all timings even of running threads are finished.
     */

    std::shared_ptr<class TraceConfig> mTraceConfig;

    class TraceData* mTraceData; // pointer to all trace data of the thread

    int mRegionId;// Reference id of region in region table.

    bool mTimeTrace;//!< set to true if timing should be done
    bool mCallTree;//!< set to true if calltree tracing should be done
    bool mVampirTrace;//!< set to true if Vampir trace should be done

    double mStartTime;//!< walltime of region start
};

class ScopedTraceRecord : private TraceRegionRecord
{
public:

    ScopedTraceRecord( const char* regionName, const char* fileName, int lno ) :

        TraceRegionRecord( regionName, fileName, lno )

    {
        SCAI_LOG_DEBUG( logger, "ScopedTraceRecord" )
        enter();
    }

    ScopedTraceRecord( const char* regionName, const int suffix_n, const char* fileName, int lno ) :

        TraceRegionRecord( regionName, suffix_n, fileName, lno )

    {
        SCAI_LOG_DEBUG( logger, "ScopedTraceRecord" )
        enter();
    }

    ~ScopedTraceRecord()
    {
        SCAI_LOG_DEBUG( logger, "~ScopedTraceRecord, call leave" )
        leave();
    }
};

} /* end namespace tracing */

} /* end namespace scai */
