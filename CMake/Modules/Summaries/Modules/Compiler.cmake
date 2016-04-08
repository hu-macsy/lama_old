heading ( "Compiler:" )

if    ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR )
    set( REQUIRED_FOUND TRUE )
else  ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR )
    set( REQUIRED_FOUND FALSE )
endif ( CXX_SUPPORTS_C11 OR SCAI_BOOST_INCLUDE_DIR )

heading2 ( "Configuration" "REQUIRED_FOUND" )
    found_message ( "C++ Compiler" "CMAKE_CXX_COMPILER" "REQUIRED" "" )
    found_message ( "with C++11 support" "CXX_SUPPORTS_C11" "REQUIRED" "" )

if    ( NOT CXX_SUPPORTS_C11 )
    emptyline()
    message ( STATUS "Either compiler supporting C++11 or Boost needed." )
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )
endif ( NOT CXX_SUPPORTS_C11 )