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

project( tracing CXX )

set ( tracing_LIBRARY tracing
      CACHE FILEPATH "trace library" )

set ( tracing_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR} 
      CACHE PATH "traing include directory" )

mark_as_advanced( tracing_LIBRARY tracing_INCLUDE_DIR )

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

add_library ( ${tracing_LIBRARY} SHARED ${CPP_SOURCES} )

include_directories( ${tracing_INCLUDE_DIR} )

### add definitions ###

find_package( Boost REQUIRED COMPONENTS )

# set ( CMAKE_INSTALL_PREFIX /home/brandes/local )

foreach ( module common logging )

   find_library( ${module}_LIBRARY 
                 NAMES ${module}
                 HINTS ${CMAKE_INSTALL_PREFIX}/lib )

   set ( ${module}_INCLUDE_DIR
            ${CMAKE_INSTALL_PREFIX}/include 
            CACHE PATH "include" )

   include_directories( ${${module}_INCLUDE_DIR} )

   target_link_libraries( ${tracing_LIBRARY} ${${module}_LIBRARY} )

endforeach()

include_directories( ${Boost_INCLUDE_DIR} )

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

add_definitions( -DLAMA_LOG_LEVEL_DEBUG )

### add library and put path to search path ###

if ( FALSE )

### VampirTrace is not easy to add for tracing, is in development

    find_package ( VampirTrace )

    message( STATUS "VAMPIRTRACE found = ${VAMPIRTRACE_FOUND}" )
    message( STATUS "MPI found = ${MPI_FOUND}" )

    if ( ${VAMPIRTRACE_FOUND} )
        message( STATUS "VAMPIRTRACE found, lib = ${VAMPIRTRACE_LIBRARIES}" )
        add_definitions( -DUSE_VAMPIRTRACE )
        include_directories( ${VAMPIRTRACE_INCLUDE_DIR} )
        target_link_libraries( ${tracing_LIBRARY} ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES} )
    endif ()

endif ()

add_subdirectory( gui )
add_subdirectory( examples )

### install ###

install ( TARGETS ${tracing_LIBRARY} DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/tracing )
