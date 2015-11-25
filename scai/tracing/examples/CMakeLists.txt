
# Explicitly enable tracing at compile time here

add_definitions( -DSCAI_TRACE_ON )

foreach ( executable TraceTest TraceCalls TraceThreads )
    add_executable( ${executable} ${executable}.cpp )
    target_link_libraries( ${executable} scai_tracing )
endforeach ( executable )

## install examples
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )
install ( DIRECTORY ./ DESTINATION share/examples/scai-tracing-${SCAI_TRACING_VERSION} FILES_MATCHING PATTERN "*.cpp" )
install ( FILES Makefile ${CMAKE_CURRENT_BINARY_DIR}/make.inc DESTINATION share/examples/scai-tracing-${SCAI_TRACING_VERSION} )