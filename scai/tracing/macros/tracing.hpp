/*
 * tracing.hpp
 *
 *  Created on: Mar 21, 2016
 *      Author: eschricker
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
