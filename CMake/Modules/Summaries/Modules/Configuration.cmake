include ( Dependencies/internal )

heading ( "Configuration Details:" )
emptyline()

set ( PROJECT_TEXT "SCAI ${PROJECT_SURNAME} Version ${SCAI_${UPPER_PROJECT_SURNAME}_VERSION}" )

if    ( ${PROJECT_NAME} MATCHES "LAMA_ALL" )
	set ( PROJECT_TEXT "${PROJECT_TEXT} ${SCAI_VERSION_NAME}" )
endif ( ${PROJECT_NAME} MATCHES "LAMA_ALL" )

indent_message ( "1" "${PROJECT_TEXT}" )
indent_message ( "1" "Build Type   : ${CMAKE_BUILD_TYPE}" )
indent_message ( "1" "Library Type : ${SCAI_LIBRARY_TYPE}" )
indent_message ( "1" "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )

list ( FIND ${UPPER_PROJECT_NAME}_INTERNAL_DEPS "scai_logging" FOUND_LOG_DEP )

if    ( ${FOUND_LOG_DEP} GREATER -1 )
	indent_message ( "1" "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" ) #opt
endif ( ${FOUND_LOG_DEP} GREATER -1 )

list ( FIND ${UPPER_PROJECT_NAME}_INTERNAL_DEPS "scai_tracing" FOUND_TRACE_DEP )
if    ( ${FOUND_TRACE_DEP} GREATER -1 )
indent_message ( "1" "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" ) #opt
endif ( ${FOUND_TRACE_DEP} GREATER -1 )

if    ( USE_CODE_COVERAGE )
    indent_message ( "1" "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )

emptyline()
