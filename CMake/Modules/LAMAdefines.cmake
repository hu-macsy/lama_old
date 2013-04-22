## LOGGING Level
#
#  Debug   : use -DLOG_LEVEL_DEBUG
#  Release : use -DLOG_LEVEL_INFO
#  
#  For serious problems: -DLOG_LEVEL_TRACE
#  For benchmarks:       -DLOG_LEVEL_OFF (or -DLOG_LEVEL_FATAL, -DLOG_LEVEL_ERROR)

if ( NOT LAMA_LOG_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_LOG_LEVEL "INFO" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    else ()
        set ( DEFAULT_LOG_LEVEL "TRACE" )
    endif ()
endif ( NOT LAMA_LOG_LEVEL )

set ( LAMA_LOG_LEVEL ${DEFAULT_LOG_LEVEL} CACHE STRING
      "Choose level of compile time logging: TRACE, DEBUG, INFO, WARN, ERROR, OFF" )

add_definitions ( -DLAMA_LOG_LEVEL_${LAMA_LOG_LEVEL} )

## ASSERT Level
#
#  Debug   : use -DASSERT_LEVEL_DEBUG
#  Release : use -DASSERT_LEVEL_ERROR
#  
#  For benchmarks:       -DASSERT_LEVEL_OFF

if ( NOT LAMA_ASSERT_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_ASSERT_LEVEL "ERROR" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    else ()
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    endif ()
endif ( NOT LAMA_ASSERT_LEVEL )

set ( LAMA_ASSERT_LEVEL ${DEFAULT_ASSERT_LEVEL} CACHE STRING
      "Choose level of ASSERT: DEBUG, ERROR, OFF" )

add_definitions ( -DLAMA_ASSERT_LEVEL_${LAMA_ASSERT_LEVEL} )

## LAMA TRACE LEVEL
#
# If TRACE is set to OFF all LAMA_REGION macros in the code are
# completely ignored. If TRACE is set to VT, regions will be traced
# (entry, exit event) for VampirTrace.

if ( NOT LAMA_TRACE_LEVEL )
    set ( DEFAULT_TRACE_LEVEL "OFF" )
endif ( NOT LAMA_TRACE_LEVEL )

set ( LAMA_TRACE_LEVEL ${DEFAULT_TRACE_LEVEL} CACHE STRING
     "Choose level of TRACE: VT (for VampirTrace), TIME(region timing), SIMPLE(simple timing) or OFF (default)" )

add_definitions( -DLAMA_TRACE_LEVEL_${LAMA_TRACE_LEVEL} )