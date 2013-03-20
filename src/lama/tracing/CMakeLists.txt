# This directory contains all relevant tracing/profiling stuff and is completely optional


### set source files ###

# Region timing is always added

set ( TRACESOURCES
        tracing/LAMABaseTracer.cpp
        tracing/TraceConfig.cpp
        tracing/RegionTable.cpp
        tracing/TraceRegionRecord.cpp
        tracing/VTInterface.cpp
    )

# for SIMPLE timing use LAMASimpleTimeTracer otherwise LAMATIMETRACER

if ( LAMA_TRACE_LEVEL STREQUAL "SIMPLE" )
    if ( Boost_THREAD_LIBRARY )
        set ( TRACESOURCES 
                tracing/LAMASimpleTimeTracer.cpp
                ${TRACESOURCES}
            )
    else ( Boost_THREAD_LIBRARY )
        message ( FATAL_ERROR "SIMPLE timing can not be used without boost thread library." )
    endif ( Boost_THREAD_LIBRARY )
endif ( LAMA_TRACE_LEVEL STREQUAL "SIMPLE" )

# PARENT_SCOPE is now only used once in this subdirectory

set ( CXXSOURCES
        ${CXXSOURCES} 
        ${TRACESOURCES} 
        PARENT_SCOPE
    )
