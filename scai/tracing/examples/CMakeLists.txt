
add_definitions( -DSCAI_TRACE_ON )


foreach ( executable TraceTest TraceCalls TraceThreads )
    add_executable( ${executable} ${executable}.cpp )
    target_link_libraries( ${executable} scai_tracing )
endforeach ( executable )
