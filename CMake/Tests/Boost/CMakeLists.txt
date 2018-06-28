###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMakeLists.txt
 # @author Thomas Brandes
 # @date 25.04.2013
###

cmake_minimum_required (VERSION 2.8 )

project( BoostExample )

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
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

