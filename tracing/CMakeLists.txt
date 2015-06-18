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
 # @brief tracing/CMakeLists.txt
 # @author Thomas Brandes
 # @date 10.06.2015
###

### Set required CMake version

cmake_minimum_required ( VERSION 2.8 )

set ( TRACING_LIBRARY scaitrace
      CACHE FILEPATH "trace library" )

set ( TRACING_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
      CACHE PATH "traing include directory" )

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

# Set default build type to RELEASE, can be changed via ccmake

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE RELEASE 
        CACHE STRING "Choose the type of build, options are: None Debug Release."
        FORCE)
endif (NOT CMAKE_BUILD_TYPE)

### Source files for tracing library

file( GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" )
file( GLOB HPP_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/tracing/*.hpp" )

### add definitions ###

# set ( CMAKE_INSTALL_PREFIX /home/brandes/local )

find_library( LOGGING_LIBRARY 
              NAMES logging 
              HINTS ${CMAKE_INSTALL_PREFIX}/lib )

set ( LOGGING_INCLUDE_DIR 
            ${CMAKE_INSTALL_PREFIX}/include 
            CACHE PATH "include" )

include_directories( ${TRACING_INCLUDE_DIR} )
include_directories( ${COMMON_INCLUDE_DIR} )
include_directories( ${LOGGING_INCLUDE_DIR} )

if ( WIN32 )
    add_definitions ( -DLOG4LAMA_COMPILING_DLL )
endif ( WIN32 )

add_definitions( -DLAMA_LOG_LEVEL_DEBUG )

### add library and put path to search path ###

add_library ( ${TRACING_LIBRARY} SHARED ${CPP_SOURCES} )
target_link_libraries( ${TRACING_LIBRARY} ${LOGGING_LIBRARY} )

find_package ( VampirTrace )
find_package ( MPI )

message( STATUS "VAMPIRTRACE found = ${VAMPIRTRACE_FOUND}" )
message( STATUS "MPI found = ${MPI_FOUND}" )

if ( ${VAMPIRTRACE_FOUND} )
    message( STATUS "VAMPIRTRACE found, lib = ${VAMPIRTRACE_LIBRARIES}" )
    add_definitions( -DUSE_VAMPIRTRACE )
    include_directories( ${VAMPIRTRACE_INCLUDE_DIR} )
    target_link_libraries( ${TRACING_LIBRARY} ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES} )
endif ()

add_subdirectory( gui )
add_subdirectory( examples )

### install ###

install ( TARGETS ${TRACING_LIBRARY} DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/tracing )
