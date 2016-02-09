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

if ( NOT DEFINED BOOST_INCLUDE_DIR )

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

find_package ( Boost ${SCAI_FIND_PACKAGE_FLAGS} COMPONENTS ${Boost_COMPONENTS} )

if    ( Boost_INCLUDE_DIR )
	set ( BOOST_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ) # for getting the module names straight
endif ( Boost_INCLUDE_DIR )

endif ( NOT DEFINED BOOST_INCLUDE_DIR )

if    ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND AND ${BUILD_TEST} )
    message ( WARNING "Boost Test Framework is missing, so BUILD_TEST is disabled!" )
    set ( BUILD_TEST FALSE )
endif ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND AND ${BUILD_TEST} )

if    ( NOT Boost_REGEX_FOUND AND ${BUILD_TEST} )
    message ( WARNING "Boost Regex is missing, so BUILD_TEST is disabled!" )
    set ( BUILD_TEST FALSE )
endif ( NOT Boost_REGEX_FOUND AND ${BUILD_TEST} )
