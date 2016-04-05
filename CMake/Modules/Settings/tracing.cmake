# SCAI_TRACING
#
# If TRACING is disabled all SCAI_REGION macros in the code are
# ignored. Otherwise performance data can be collected
# where configuration is set at runtime via SCAI_TRACE.

include ( Functions/checkValue )

set ( SCAI_TRACING ON )
checkValue ( ${SCAI_TRACING} "${TRUE_FALSE_CHOICES}" )

if ( SCAI_TRACING )
    set ( SCAI_TRACING_FLAG "SCAI_TRACE_ON" )
else ( SCAI_TRACING )
    set ( SCAI_TRACING_FLAG "SCAI_TRACE_OFF" )
endif ( SCAI_TRACING )
