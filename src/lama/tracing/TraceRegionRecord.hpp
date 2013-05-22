/**
 * @file TraceRegionRecord.hpp
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
 * @brief Definition of class for tracing that measures inclusive and exclusive
 *        walltimes for regions defined by directive LAMA_REGION.
 * @author Thomas Brandes
 * @date 01.09.2011
 * @since 1.0.0
 */
#ifndef LAMA_TRACING_TRACE_REGION_RECORD_HPP_
#define LAMA_TRACING_TRACE_REGION_RECORD_HPP_

// for dll_import
#include <lama/config.hpp>

// base classes
#include <lama/tracing/LAMABaseTracer.hpp>

// logging
#include <logging/logging.hpp>

// boost
#include <boost/shared_ptr.hpp>

/** This class is a helper class for tracing a scope in C++.
 *
 *  The constructor of an object creates a start entry and the destructor the stop entry
 *  for trace files. The use of this class is much safer than using calls of routines start and
 *  stop explicitly as it works also very well in case of exception.
 */

namespace tracing
{

class TraceRegionRecord: public LAMABaseTracer
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

private:

    LAMA_LOG_DECL_STATIC_LOGGER( logger )

    /** Common routine used in all constructors. */

    void enter( const char* regionName, const char* file, int lno );

    /** Each region timing keeps a shared pointer to the configuration.
     *  By this way it is guaranteed that timer information is only printed
     *  when all timings even of running threads are finished.
     */

    boost::shared_ptr<class TraceConfig> mTraceConfig;

    class RegionTable* mRegionTable; // pointer to thread region time table

    int mRegionId;// Reference id of region in region table.

    bool mTimeTrace;//!< set to true if timing should be done
    bool mVampirTrace;//!< set to true if Vampir trace should be done

    double mStartTime;//!< walltime of region start
};

} // namespace

#endif // LAMA_TRACING_TRACE_REGION_RECORD_HPP_
