/**
 * @file tracing.hpp
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
 * @brief Definition of macros for tracing/profiling
 * @author Lauretta Schubert, Thomas Brandes
 * @date 12.08.2011
 * @since 1.0.0
 */

#pragma once

#if defined( LAMA_TRACE_ON )

#include <tracing/TraceRegionRecord.hpp>
#include <tracing/TraceConfig.hpp>

#define LAMA_REGION( name ) tracing::ScopedTraceRecord LAMA_Trc__( name, __FILE__, __LINE__ );
#define LAMA_REGION_N( name, n ) tracing::ScopedTraceRecord LAMA_Trc__( name, n, __FILE__, __LINE__ );
#define LAMA_REGION_START( name ) tracing::TraceRegionRecord::start( name, __FILE__, __LINE__ );
#define LAMA_REGION_END( name ) tracing::TraceRegionRecord::stop( name );
#define LAMA_TRACE_SCOPE( flag ) tracing::TraceConfig::TraceScope LAMA_Scp__( flag );
#define LAMA_TIMETRACER( name ) tracing::TraceRegionRecord::spentLast( name );

#elif defined( LAMA_TRACE_OFF )

#define LAMA_REGION( name )
#define LAMA_REGION_START( name )
#define LAMA_REGION_END( name )
#define LAMA_REGION_N( name, n )
#define LAMA_TRACE_SCOPE( flag )
#define LAMA_TIMETRACER( name ) 0.0

#else

//Macro LAMA_REGION should also be defined in case of error for convenience with Eclipse.
#define LAMA_REGION( name )
#define LAMA_REGION_START( name )
#define LAMA_REGION_END( name )
#define LAMA_REGION_N( name, n )
#define LAMA_TRACE_SCOPE( flag )
#define LAMA_TIMETRACER( name ) 0.0
#error "Must define LAMA_TRACE_xxx with xxx = ON, or OFF"

#endif
