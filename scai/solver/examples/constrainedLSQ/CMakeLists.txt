###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief CMake configuration file for tutorial examples
 # @author Thomas Brandes
 # @date 25.07.2017
###

project ( LeastSquare CXX )

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

add_executable ( prepare prepare.cpp )
target_link_libraries( prepare ${SCAI_LIBRARIES} )

add_executable ( lsbc lsbc.cpp ConstrainedLeastSquares.cpp )
target_link_libraries( lsbc ${SCAI_LIBRARIES} )

add_executable ( makeImage makeImage.cpp )
target_link_libraries( makeImage ${SCAI_LIBRARIES} )

add_executable ( doAll doAll.cpp ConstrainedLeastSquares.cpp )
target_link_libraries( doAll ${SCAI_LIBRARIES} )

