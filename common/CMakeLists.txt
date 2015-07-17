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
 # @brief CMake configuration file for common library
 # @author Thomas Brandes
 # @date 09.06.2015
###

### Set required CMake version

cmake_minimum_required ( VERSION 2.8 )

project ( common CXX )

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

# Set default build type to RELEASE, can be changed via ccmake

if ( NOT CMAKE_BUILD_TYPE )
   set ( CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif ( NOT CMAKE_BUILD_TYPE )

include ( SetCPPFlags )

### Define variables needed for common library 
### Variables are put in the cache so they might be used in other scopes

set ( common_LIBRARY common 
      CACHE FILEPATH "common library" )

set ( common_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
      CACHE PATH "common include directory" )

mark_as_advanced ( common_INCLUDE_DIR common_LIBRARY )

###  Here we use PThread library for threads
###  Note: FindThreads in CMake is available as Module, but is buggy, needs update of CheckIncludeFiles.cmake

find_library ( PTHREADS_LIBRARY NAMES pthread pthreads )

file ( GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" )
file ( GLOB HPP_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/common/*.hpp" )

include_directories ( ${common_INCLUDE_DIR} )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### add library and put path to search path ###

add_library ( ${common_LIBRARY} SHARED ${CPP_SOURCES} )
target_link_libraries( ${common_LIBRARY} ${PTHREADS_LIBRARY} )

add_subdirectory ( examples )

### install ###

install ( TARGETS ${common_LIBRARY} DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/common )
