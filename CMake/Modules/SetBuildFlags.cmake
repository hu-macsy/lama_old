###
 # @file Variables.cmake
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
 # @brief Important CMake variable definitions
 # @author Jan Ecker
 # @date 16.04.2013
 # @since 1.0.0
###

include ( Functions/checkValue )

# Check if verbose mode for CMAKE is selected
if    ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    set ( SCAI_FIND_PACKAGE_FLAGS )
else  ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    set ( SCAI_FIND_PACKAGE_FLAGS QUIET )
endif ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )

## BUILDTYPE

# Set default build type, will be used for all projects
# Note: can be changed at any time via CCMAKE

# Choose Default CMAKE_BUILD_TYPE
if ( NOT CMAKE_BUILD_TYPE )
	set ( CMAKE_BUILD_TYPE_OPTIONS None Debug Release RelWithDebInfo MinSizeRel ) 
    # Can be: (RelWithDebInfo)
    set ( CMAKE_BUILD_TYPE Debug CACHE STRING 
        "Choose the type of build, options are: ${CMAKE_BUILD_TYPE_OPTIONS}." FORCE )
	checkValue ( ${CMAKE_BUILD_TYPE} "${CMAKE_BUILD_TYPE_OPTIONS}" )
    message ( STATUS "Build type is set to " ${CMAKE_BUILD_TYPE} )
endif ( NOT CMAKE_BUILD_TYPE )

# Choose Doc type
if ( NOT SCAI_DOC_TYPE )
    set ( SCAI_DOC_TYPE_OPTIONS html json ) 
    set ( SCAI_DOC_TYPE html CACHE STRING 
        "Choose the type of documentation, options are: ${SCAI_DOC_TYPE_OPTIONS}." FORCE )
    checkValue ( ${SCAI_DOC_TYPE} "${SCAI_DOC_TYPE_OPTIONS}" )
    message ( STATUS "Doc type is set to " ${CMAKE_BUILD_TYPE} )

    if    ( SCAI_DOC_TYPE STREQUAL json )
        set ( DOC_EXTENTSION "fjson" )
    else  ( SCAI_DOC_TYPE STREQUAL json )
        set ( DOC_EXTENTSION "html" )
    endif ( SCAI_DOC_TYPE STREQUAL json )
endif ( NOT SCAI_DOC_TYPE )

## Check if lama should be build static or shared

# default: build shared library
set ( SCAI_LIBRARY_TYPE_OPTIONS STATIC SHARED )
if    ( NOT SCAI_LIBRARY_TYPE )
	set ( SCAI_LIBRARY_TYPE SHARED CACHE STRING "Choose the type of linking: ${SCAI_LIBRARY_TYPE_OPTIONS}" )
else  ( NOT SCAI_LIBRARY_TYPE ) 	
    set ( SCAI_LIBRARY_TYPE STATIC CACHE STRING "Choose the type of linking: ${SCAI_LIBRARY_TYPE_OPTIONS}" )
endif ( NOT SCAI_LIBRARY_TYPE )
checkValue ( ${SCAI_LIBRARY_TYPE} "${SCAI_LIBRARY_TYPE_OPTIONS}" )


set ( TRUE_FALSE_CHOICE ON OFF )

# Set cache variable
set ( BUILD_TEST TRUE CACHE BOOL "Enable / Disable building of tests" )
checkValue ( ${BUILD_TEST} "${TRUE_FALSE_CHOICE}" )

## Check if lama should be build for code coverage
set ( USE_CODE_COVERAGE FALSE CACHE BOOL "Enable / Disable use of Code Coverage" )
checkValue ( ${USE_CODE_COVERAGE} "${TRUE_FALSE_CHOICE}" )
