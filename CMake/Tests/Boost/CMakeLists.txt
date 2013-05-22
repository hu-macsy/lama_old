###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2013
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
 # @brief CMakeLists.txt
 # @author Thomas Brandes
 # @date 25.04.2013
 # $Id$
###

cmake_minimum_required (VERSION 2.8 )

project( BoostExample )

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set ( Boost_COMPONENTS thread )

set ( Boost_DEBUG TRUE )
set ( Boost_DETAILED_FAILURE_MSG TRUE )

# Alternative on WIN32

set ( Boost_USE_STATIC_LIBS ON )

# set ( Boost_USE_MULTITHREADED ON )
# set ( Boost_USE_STATIC_RUNTIME OFF )

# Note: set (environment) variable BOOST_ROOT to find installed Boost

find_package ( Boost COMPONENTS ${Boost_COMPONENTS} )

if ( NOT Boost_THREAD_LIBRARY )
    message ( FATAL_ERROR "Boost thread library not found" )
endif ( NOT Boost_THREAD_LIBRARY )

if ( Boost_FOUND )
    message( STATUS "Boost found" )
    message( STATUS "Boost_THREAD_LIBRARY=${Boost_THREAD_LIBRARY}" )
    message( STATUS "Boost_LIBRARY_DIRS=${Boost_LIBRARY_DIRS}" )
    include_directories( ${Boost_INCLUDE_DIR} )

    # For static libraries on Windows it is likely that others must be added
    # (e.g. data_time, system); so add library directory
    # must appear before executable is defined

    link_directories( ${Boost_LIBRARY_DIRS} )

    add_executable( example example.cpp )

    target_link_libraries( example ${Boost_THREAD_LIBRARY} )

else( Boost_FOUND )

    message( STATUS "Boost not found" )

endif( Boost_FOUND )

