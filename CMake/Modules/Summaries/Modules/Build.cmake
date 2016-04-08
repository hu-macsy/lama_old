heading ( "Build options:" "" )

# EXAMPLES
heading3 ( "Examples" "BUILD_EXAMPLES" )

# TEST
if    ( NOT ( ( ${PROJECT_NAME} MATCHES "scai_logging" ) OR ( ${PROJECT_NAME} MATCHES "scai_tracing" ) ) )

	heading3 ( "Test" "BOOST_TEST_ENABLED" )
	    found_message ( "Boost Unit Test" "Boost_UNIT_TEST_FRAMEWORK_FOUND" "OPTIONAL" "Version ${BOOST_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
	    found_message ( "Boost Regex"     "Boost_REGEX_FOUND"               "OPTIONAL" "Version ${BOOST_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

endif ( NOT ( ( ${PROJECT_NAME} MATCHES "scai_logging" ) OR ( ${PROJECT_NAME} MATCHES "scai_tracing" ) ) )

# DOC
heading3 ( "Documentation" "DOC_ENABLED" )
    found_message ( "Sphinx" "SPHINX_FOUND" "OPTIONAL" "Version ${Sphinx_VERSION_STRING} with ${Sphinx-build_EXECUTABLE}" )
    found_message ( "Doxygen" "DOXYGEN_FOUND" "OPTIONAL" "Version ${DOXYGEN_VERSION} with ${DOXYGEN_EXECUTABLE}" )