/**
 * @file tracing.hpp
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
 * @brief Definition of macros for tracing/profiling
 * @author Lauretta Schubert, Thomas Brandes
 * @date 12.08.2011
 * $Id$
 */
#ifndef LAMA_TRACING_HPP_
#define LAMA_TRACING_HPP_

/* ToDo: use it from CUDA directory
 #ifdef LAMA_BUILD_CUDA
 #include <lama/tracing/CUDATracerHelper.hpp>
 #endif
 */

#if defined( LAMA_TRACE_LEVEL_VT ) || defined( LAMA_TRACE_LEVEL_TIME )

#include <lama/tracing/TraceRegionRecord.hpp>

#define LAMA_REGION( name ) tracing::TraceRegionRecord LAMA_Trc__( name, __FILE__, __LINE__ )

#define LAMA_REGION_N( name, n ) tracing::TraceRegionRecord LAMA_Trc__( name, n, __FILE__, __LINE__ )

#define LAMA_REGION_START( name ) tracing::TraceRegionRecord::start( name, __FILE__, __LINE__ )

#define LAMA_REGION_END( name ) tracing::TraceRegionRecord::stop( name )

/* ToDo: use it from CUDA subdirectory
 #ifdef LAMA_BUILD_CUDA
 #define LAMA_CUDAREGION( name, syncToken ) \
        CUDATracerHelper<tracing::TraceRegionRecord> LAMA_Trc__( name, __FILE__, __LINE__, syncToken )
 #else
 #define LAMA_CUDAREGION( name, syncToken )
 #endif

 #ifdef LAMA_BUILD_OPENCL
 //TODO: OpenCLTracerHelper
 #define LAMA_OPENCLREGION( name, syncToken )
 #else
 #define LAMA_OPENCLREGION( name, syncToken )
 #endif
 */

#define LAMA_TIMETRACER( name ) tracing::TraceRegionRecord::spentLast( name )

#elif defined( LAMA_TRACE_LEVEL_SIMPLE )

#include <lama/tracing/LAMASimpleTimeTracer.hpp>

#define LAMA_REGION( name ) LAMASimpleTimeTracer LAMA_Trc__( name, __FILE__, __LINE__ )

#define LAMA_REGION_N( name, n )

#define LAMA_REGION_START( name )

#define LAMA_REGION_END( name )

#ifdef LAMA_BUILD_CUDA
#define LAMA_CUDAREGION( name, syncToken ) \
    CUDATracerHelper<LAMASimpleTimeTracer> LAMA_Trc__( name, __FILE__, __LINE__, syncToken )
#else
#define LAMA_CUDAREGION( name, syncToken )
#endif

#ifdef LAMA_BUILD_OPENCL
//TODO: OpenCLTracerHelper
#define LAMA_OPENCLREGION( name, syncToken )
#else
#define LAMA_OPENCLREGION( name, syncToken )
#endif

#define LAMA_TIMETRACER( name ) LAMASimpleTimeTracer::spentLast( name )

#elif defined( LAMA_TRACE_LEVEL_OFF )

#define LAMA_REGION( name )

#define LAMA_REGION_START( name )

#define LAMA_REGION_END( name )

#define LAMA_REGION_N( name, n )

#define LAMA_CUDAREGION( name, syncToken )

#define LAMA_OPENCLREGION( name, syncToken )

#define LAMA_TIMETRACER( name ) 0.0

#else

//Macro LAMA_REGION should also be defined in case of error for convenience with Eclipse.

#define LAMA_REGION( name )

#define LAMA_REGION_START( name )

#define LAMA_REGION_END( name )

#define LAMA_REGION_N( name, n )

#define LAMA_CUDAREGION( name, syncToken )

#define LAMA_OPENCLREGION( name, syncToken )

#define LAMA_TIMETRACER( name ) 0.0

#error "Must define LAMA_TRACE_LEVEL_xxx with xxx = VT, TIME, SIMPLE, or OFF"

#endif

#endif // LAMA_TRACING_HPP_
