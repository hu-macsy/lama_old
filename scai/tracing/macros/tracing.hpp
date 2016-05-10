/**
 * @file macros/tracing.hpp
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
 * @brief Definition of the tracing macros.
 * @author Lauretta Schubert, Thomas Brandes
 * @date 12.08.2011
 */

#pragma once

#if defined( SCAI_TRACE_ON )

    #include <scai/tracing/TraceRegionRecord.hpp>
    #include <scai/tracing/TraceConfig.hpp>

    #define SCAI_REGION( name ) scai::tracing::ScopedTraceRecord SCAI_Trc__( name, __FILE__, __LINE__ );
    #define SCAI_REGION_N( name, n ) scai::tracing::ScopedTraceRecord SCAI_Trc__( name, n, __FILE__, __LINE__ );
    #define SCAI_REGION_START( name ) scai::tracing::TraceRegionRecord::start( name, __FILE__, __LINE__ );
    #define SCAI_REGION_END( name ) scai::tracing::TraceRegionRecord::stop( name );
    #define SCAI_TRACE_SCOPE( flag ) scai::tracing::TraceConfig::TraceScope SCAI_Scp__( flag );

#elif defined( SCAI_TRACE_OFF )

    #define SCAI_REGION( name )
    #define SCAI_REGION_START( name )
    #define SCAI_REGION_END( name )
    #define SCAI_REGION_N( name, n )
    #define SCAI_TRACE_SCOPE( flag )

#else

    //Macro SCAI_REGION should also be defined in case of error for convenience with Eclipse.
    #define SCAI_REGION( name )
    #define SCAI_REGION_START( name )
    #define SCAI_REGION_END( name )
    #define SCAI_REGION_N( name, n )
    #define SCAI_TRACE_SCOPE( flag )
    #error "Must define SCAI_TRACE_xxx with xxx = ON, or OFF"

#endif
