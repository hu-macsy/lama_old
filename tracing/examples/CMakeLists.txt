
add_definitions( -DLAMA_TRACE_LEVEL_TIME )

add_executable( TraceTest TraceTest.cpp )
target_link_libraries( TraceTest ${TRACING_LIBRARY} )
