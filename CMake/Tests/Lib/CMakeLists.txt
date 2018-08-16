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

################################################
#   CMake example for simple C++ program       #
################################################

cmake_minimum_required( VERSION 2.8 )

## PROJECT NAME and supported LANGUAGES

project ( LibTest CXX )

################################################
#   Simple C++ program                         #
################################################

# default is building shared libraries, might be changed

set ( BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries instead of static ones." )

if ( BUILD_SHARED_LIBS )

    set ( LIBRARY_MODE SHARED )

    add_definitions( -DDLL_LINKING )

else ()

    set ( LIBRARY_MODE STATIC )

endif ()

add_library( rectangle ${LIBRARY_MODE} Rectangle.cpp Rectangle.hpp )

add_executable( example example.cpp )

target_link_libraries( example rectangle )
