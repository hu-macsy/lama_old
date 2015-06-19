
add_definitions( -DLAMA_TRACE_ON )

add_executable( TraceTest TraceTest.cpp )
target_link_libraries( TraceTest ${TRACING_LIBRARY} )

add_executable( TraceCalls TraceCalls.cpp )
target_link_libraries( TraceCalls ${TRACING_LIBRARY} )
