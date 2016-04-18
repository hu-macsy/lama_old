# SCAI_TRACING
#
# If TRACING is disabled all SCAI_REGION macros in the code are
# ignored. Otherwise performance data can be collected
# where configuration is set at runtime via SCAI_TRACE.

include ( Functions/checkValue )

if    ( CMAKE_BUILD_TYPE STREQUAL "Release" )
	set ( SCAI_TRACING OFF )
else  ( CMAKE_BUILD_TYPE STREQUAL "Release" )
	set ( SCAI_TRACING ON )
endif ( CMAKE_BUILD_TYPE STREQUAL "Release" )
checkValue ( ${SCAI_TRACING} "${TRUE_FALSE_CHOICES}" )
set ( SCAI_TRACING ${SCAI_TRACING} CACHE BOOL "Enable / Disable tracing of regions for performance analysis" )

if ( SCAI_TRACING )
    set ( SCAI_TRACING_FLAG "SCAI_TRACE_ON" )
else ( SCAI_TRACING )
    set ( SCAI_TRACING_FLAG "SCAI_TRACE_OFF" )
endif ( SCAI_TRACING )
