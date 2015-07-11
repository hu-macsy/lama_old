
add_definitions( -DLAMA_TRACE_ON )


foreach ( executable TraceTest TraceCalls TraceThreads )
    add_executable( ${executable} ${executable}.cpp )
    target_link_libraries( ${executable} ${tracing_LIBRARY} )
endforeach ( executable )
