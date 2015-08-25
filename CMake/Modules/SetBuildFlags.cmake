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

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

get_property ( FIND_LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS )
message ( STATUS "FindLib64: " ${FIND_LIB64} )

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )

# Define default library type as SHARED
if ( NOT DEFINED BUILD_SHARED_LIBS )
    set ( BUILD_SHARED_LIBS TRUE )
endif( NOT DEFINED BUILD_SHARED_LIBS )

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
endif ( NOT CMAKE_BUILD_TYPE )

message ( STATUS "Build type is set to " ${CMAKE_BUILD_TYPE} )

## Check if lama should be build static or shared
## Check if cache variable is already set
if    ( DEFINED SCAI_BUILD_SHARED )
    set ( SCAI_BUILD_SHARED ${SCAI_BUILD_SHARED} CACHE BOOL "Enable / Disable dynamic linking of libraries" )
else  ( DEFINED SCAI_BUILD_SHARED )
    # Set cache variable
    set ( SCAI_BUILD_SHARED TRUE CACHE BOOL "Enable / Disable dynamic linking of libraries" )
endif ( DEFINED SCAI_BUILD_SHARED )

if    ( SCAI_BUILD_SHARED )
    set ( SCAI_LIBRARY_TYPE SHARED )
else  ( SCAI_BUILD_SHARED )
    set ( SCAI_LIBRARY_TYPE STATIC )
endif ( SCAI_BUILD_SHARED )

set ( USE_CODE_COVERAGE FALSE CACHE BOOL "Enable / Disable use of Code Coverage" )