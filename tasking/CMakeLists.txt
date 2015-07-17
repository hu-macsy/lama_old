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
 # @brief tasking/CMakeLists.txt
 # @author Thomas Brandes
 # @date 10.06.2015
###

### Set required CMake version

cmake_minimum_required ( VERSION 2.8 )

project( tasking CXX )

set ( tasking_LIBRARY tasking
      CACHE FILEPATH "tasking library" )

set ( tasking_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
      CACHE PATH "tasking include directory" )

mark_as_advanced( tasking_LIBRARY tasking_INCLUDE_DIR )

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

# Set default build type to RELEASE, can be changed via ccmake

if ( NOT CMAKE_BUILD_TYPE )
   set ( CMAKE_BUILD_TYPE RELEASE 
        CACHE STRING "Choose the type of build, options are: None Debug Release."
        FORCE)
endif ( NOT CMAKE_BUILD_TYPE )

include ( SetCPPFlags )

### Source files for tasking library

file( GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp" )
file( GLOB HPP_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/tasking/*.hpp" )

### add library and put path to search path ###

add_library ( ${tasking_LIBRARY} SHARED ${CPP_SOURCES} )

include_directories( ${tasking_INCLUDE_DIR} )

### add definitions ###

# set ( CMAKE_INSTALL_PREFIX /home/brandes/local )

find_package( Boost REQUIRED COMPONENTS )

foreach ( module common logging tracing )

   find_library( ${module}_LIBRARY 
                 NAMES ${module}
                 HINTS ${CMAKE_INSTALL_PREFIX}/lib )

   set ( ${module}_INCLUDE_DIR
            ${CMAKE_INSTALL_PREFIX}/include 
            CACHE PATH "include" )

   include_directories( ${${module}_INCLUDE_DIR} )

   target_link_libraries( ${tasking_LIBRARY} ${${module}_LIBRARY} )

endforeach()

include_directories( ${Boost_INCLUDE_DIR} )

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

add_definitions( -DLAMA_LOG_LEVEL_DEBUG )
add_definitions( -DLAMA_TRACE_ON )

add_subdirectory( examples )

### install ###

install ( TARGETS ${tasking_LIBRARY} DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/tasking )
