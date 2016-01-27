###
 # @file PackageBoost.cmake
 #
 # @license
 # Copyright (c) 2009-2013
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief findPackage and configuration of Boost
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

### BOOST_<lib>_FOUND    - if Boost component is found
### BOOST_INCLUDE_DIR    - Boost include directory
### Boost_<lib>_LIBRARY  - Boost component library
### SCAI_BOOST_LIBRARIES - all found BOOST libraries out of the searched component

### Boost_USE_STATIC_LIBS  ( default is OFF )
### 
###    - should be set to OFF, default ( will use dynamic Boost libraries if available )
###    - should only be set to ON  if only static BOOST libraries are availabe 
###    - On Linux systems static 64-bit libraries must have been compiled with -fPIC

# set ( Boost_USE_STATIC_LIBS OFF )
# set ( Boost_USE_MULTITHREADED OFF )

if ( NOT DEFINED Boost_INCLUDE_DIR )

if    ( WIN32 )
    message ( STATUS "Setting special Boost options on Windows" )
    #set ( Boost_USE_STATIC_LIBS ON )
    set ( Boost_USE_MULTITHREADED ON )
endif ( WIN32 )

# Finds packages with custom search options 

set ( Boost_COMPONENTS unit_test_framework regex thread system )

# FindBoost Debug options comment
if    ( SCAI_CMAKE_VERBOSE )
    set ( Boost_DEBUG TRUE )
    set ( Boost_DETAILED_FAILURE_MSG TRUE )
endif ( SCAI_CMAKE_VERBOSE )

# Find Boost 

#set ( SCAI_CMAKE_VERBOSE QUIET )

find_package ( Boost ${SCAI_FIND_PACKAGE_FLAGS} COMPONENTS ${Boost_COMPONENTS} )

if    ( Boost_INCLUDE_DIR )
	set ( BOOST_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ) # for getting the module names straight
endif ( Boost_INCLUDE_DIR )

# Note: we use BOOST_INCLUDE_DIR, Boost_<lib>_FOUND, Boost_<lib>_LIBRARY, but
#       not Boost_FOUND, as it is false if some optional libraries are missing

# Boost: include directory is mandatory ( LAMA uses shared pointer, function )

if    ( BOOST_INCLUDE_DIR )
    # Boost_PATH should be same as BOOST_ROOT
    get_filename_component ( Boost_PATH "${BOOST_INCLUDE_DIR}" PATH )
else  ( BOOST_INCLUDE_DIR )
    message ( STATUS "Boost (include directory) not found: give hint by environment variable BOOST_ROOT" ) 
endif ( BOOST_INCLUDE_DIR )

# check status of each Boost component and conclude all needed Boost libraries
set ( SCAI_BOOST_LIBRARIES "" )

foreach    ( lib ${Boost_COMPONENTS} )
    string ( TOUPPER ${lib} libname )
    set ( libname "Boost_${libname}_LIBRARY" )
    # libname: variable that contains the library for the boost component
    # TODO: SUMMARY   
    #message ( STATUS "${libname} = ${${libname}}" )
    if    ( ${libname} )
        # library found, make sure it belongs to same version of Boost
        if    ( "${${libname}}" MATCHES "${Boost_PATH}*" )
        	list ( APPEND SCAI_BOOST_LIBRARIES ${${libname}} )
            # do nothing
        else  ( "${${libname}}" MATCHES "${Boost_PATH}*" )
            message ( FATAL_ERROR "${${libname}} illegal, not in ${Boost_PATH}" )
        endif ( "${${libname}}" MATCHES "${Boost_PATH}*" )
    endif ( ${libname} )
endforeach ( lib ${Boost_COMPONENTS} )

# Check boost versions
# TODO: RECHECK
if    ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )
    set ( FOUND_BOOST_TEST TRUE )
else  ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )
    if    ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND )
        message ( WARNING "Not building tests because Boost unit test framework is missing." )
    endif ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND )
    
    if    ( NOT Boost_REGEX_FOUND )
        message ( WARNING "Not building tests because Boost regex is missing." )
    endif ( NOT Boost_REGEX_FOUND )
    
    if ( ${Boost_VERSION} LESS "104100" )
        message ( WARNING "Not building tests because Boost is to old: ${Boost_VERSION}." )
    endif ( ${Boost_VERSION} LESS "104100" )
endif ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )


# if use of package is enabled
if    ( ${BUILD_TEST} AND NOT ${FOUND_BOOST_TEST} )
    # if package is enabled, but not found: ERROR!
    message ( STATUS "Boost Test Framework or Bost Regex missing, but tests are enabled!" )
endif ( ${BUILD_TEST} AND NOT ${FOUND_BOOST_TEST} )

endif ( NOT DEFINED Boost_INCLUDE_DIR )
