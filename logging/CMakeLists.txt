###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @brief CMake configuration file for Logging library
 # @author Thomas Brandes
 # @date 09.06.2015
 # @since 2.0.0
###

### Set required CMake version

cmake_minimum_required ( VERSION 2.8 )

project( logging CXX )

### Define variables needed for logging library 
### Variables are put in the cache so they might be used in other scopes

set ( LOGGING_LIBRARY logging 
      CACHE FILEPATH "logging library" )

set ( LOGGING_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
      CACHE PATH "logging include directory" )

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

# Set default build type to RELEASE, can be changed via ccmake

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

### Boost thread required

set ( Boost_DEBUG TRUE )
set ( Boost_DETAILED_FAILURE_MSG TRUE )

find_package ( Boost COMPONENTS thread REQUIRED )

# Common library of other project is used here
# It might be either installed or variables are already in CACHE

find_library( COMMON_LIBRARY 
              NAMES common
              HINTS ${CMAKE_INSTALL_PREFIX}/lib )

set ( COMMON_INCLUDE_DIR ${CMAKE_INSTALL_PREFIX}/include CACHE PATH "include" )

### Source files for logging library

file( GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" )
file( GLOB HPP_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/logging/*.hpp" )

### add definitions ###

# GenLogger uses also logging

add_definitions ( -DLAMA_LOG_LEVEL_DEBUG )

include_directories( ${Boost_INCLUDE_DIR} )
include_directories( ${COMMON_INCLUDE_DIR} )
include_directories( ${LOGGING_INCLUDE_DIR} )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### add library and put path to search path ###

add_library ( ${LOGGING_LIBRARY} SHARED ${CPP_SOURCES} )
target_link_libraries( ${LOGGING_LIBRARY} ${COMMON_LIBRARY} ${Boost_THREAD_LIBRARY} )

add_subdirectory( examples )

### install ###

install ( TARGETS ${LOGGING_LIBRARY} DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/logging )
