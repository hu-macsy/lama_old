# set ( Boost_USE_MULTITHREADED OFF )

if ( WIN32 )
    message ( STATUS "Setting special Boost options on Windows" )
    #set ( Boost_USE_STATIC_LIBS ON )
    set ( Boost_USE_MULTITHREADED ON )
    set ( Boost_USE_STATIC_RUNTIME OFF )
    #add_definitions ( -DBOOST_ALL_NO_LIB )
endif ( WIN32 )

# Finds packages with custom search options 

set ( Boost_COMPONENTS thread unit_test_framework regex )

# FindBoost Debug options comment
if ( LAMA_DEBUG_CMAKE )
    set ( Boost_DEBUG TRUE )
    set ( Boost_DETAILED_FAILURE_MSG TRUE )
endif( LAMA_DEBUG_CMAKE )

# Find Boost 

find_package ( Boost COMPONENTS ${Boost_COMPONENTS} QUIET )

# Note: we use Boost_INCLUDE_DIR, Boost_<lib>_FOUND, Boost_<lib>_LIBRARY, but
#       not Boost_FOUND, as it is false if some optional libraries are missing

# Boost: include directory is mandatory ( LAMA uses shared pointer, function )

if ( Boost_INCLUDE_DIR )
# TODO: SUMMARY
    #message ( STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}" )
    get_filename_component ( Boost_PATH ${Boost_INCLUDE_DIR} PATH )
# TODO: SUMMARY
    #message ( STATUS "Boost_PATH = ${Boost_PATH}" )
    # Boost_PATH should be same as BOOST_ROOT
else ( Boost_INCLUDE_DIR )
    message ( FATAL_ERROR "Boost (include directory) not found: give hint by environment variable BOOST_ROOT" ) 
endif ( Boost_INCLUDE_DIR )

# check status of each Boost component

foreach ( lib ${Boost_COMPONENTS} )
   string ( TOUPPER ${lib} libname )
   set ( libname "Boost_${libname}_LIBRARY" )
   # libname: variable that contains the library for the boost component
# TODO: SUMMARY   
   #message ( STATUS "${libname} = ${${libname}}" )
   if ( ${libname} )
       # library found, make sure it belongs to same version of Boost
       if ( "${${libname}}" MATCHES "${Boost_PATH}*" )
           #
       else ( "${${libname}}" MATCHES "${Boost_PATH}*" )
           message ( FATAL_ERROR "${${libname}} illegal, not in ${Boost_PATH}" )
       endif ( "${${libname}}" MATCHES "${Boost_PATH}*" )
   endif ( ${libname} )
endforeach ( lib ${Boost_COMPONENTS} )

if ( NOT Boost_THREAD_LIBRARY )
    message ( FATAL_ERROR "Boost thread library not found" )
endif ( NOT Boost_THREAD_LIBRARY )


# Check boost versions
# TODO: RECHECK
if ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )
    set ( LAMA_BUILD_TEST TRUE )
else ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )
    if ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND )
       message ( WARNING "Not building tests because Boost unit test framework is missing." )
    endif ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND )
    if ( NOT Boost_REGEX_FOUND )
       message ( WARNING "Not building tests because Boost regex is missing." )
    endif ( NOT Boost_REGEX_FOUND )
    if ( ${Boost_VERSION} LESS "104100" )
       message ( WARNING "Not building tests because Boost is to old: ${Boost_VERSION}." )
    endif ( ${Boost_VERSION} LESS "104100" )
endif ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )



