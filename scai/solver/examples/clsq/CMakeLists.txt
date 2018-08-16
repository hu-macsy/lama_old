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
 # @brief CMake configuration file for tutorial examples
 # @author Thomas Brandes
 # @date 25.07.2017
###

project ( LeastSquare )

cmake_minimum_required ( VERSION 2.8 )

find_package( SCAI REQUIRED PATHS $ENV{SCAI_ROOT} )

message ( STATUS "SCAI_DEFINITIONS=${SCAI_DEFINITIONS}" )
message ( STATUS "SCAI_INCLUDE_DIRS=${SCAI_INCLUDE_DIRS}" )
message ( STATUS "SCAI_CXX_FLAGS=${SCAI_CXX_FLAGS}" )
message ( STATUS "SCAI_LIBRARIES=${SCAI_LIBRARIES}" )

add_definitions ( ${SCAI_DEFINITIONS} )

include_directories( ${SCAI_INCLUDE_DIRS} )

## Take over OpenMP and C++11 flags 

set ( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${SCAI_CXX_FLAGS} )

add_executable ( lsbc lsbc.cpp ConstrainedLeastSquares.cpp )
target_link_libraries( lsbc ${SCAI_LIBRARIES} )
