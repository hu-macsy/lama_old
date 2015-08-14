/**
 * @file TraceRegionRecord.hpp
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
 * @brief Definition of class for tracing that measures inclusive and exclusive
 *        walltimes for regions defined by directive LAMA_REGION.
 * @author Thomas Brandes
 * @date 11.06.2015
 */

#pragma once

// for dll_import
#include <scai/common/config.hpp>
#include <scai/logging.hpp>

#include <scai/common/shared_ptr.hpp>

/** Namespace for all data structures used in tracing library. */

namespace tracing
{

/** This class is a helper class for tracing a scope in C++.
 *
 *  The constructor of an object creates a start entry and the destructor the stop entry
 *  for trace files. The use of this class is much safer than using calls of routines start and
 *  stop explicitly as it works also very well in case of exception.
 */

class TraceRegionRecord
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

    common::shared_ptr<class TraceConfig> mTraceConfig;

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

}  // namescape

