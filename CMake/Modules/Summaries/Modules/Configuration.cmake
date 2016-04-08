
heading ( "Configuration Details:" )

indent_message ( "1" "LAMA (ALL) Version : ${SCAI_LAMA_ALL_VERSION} ${SCAI_VERSION_NAME}" )
indent_message ( "1" "Build Type   : ${CMAKE_BUILD_TYPE}" )
indent_message ( "1" "Library Type : ${SCAI_LIBRARY_TYPE}" )
indent_message ( "1" "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
indent_message ( "1" "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" ) #opt
indent_message ( "1" "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" ) #opt
if    ( USE_CODE_COVERAGE )
    indent_message ( "1" "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )
