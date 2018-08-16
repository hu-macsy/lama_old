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
 # @date 18.04.2018
###

################################################
#   CMake example for simple C++ program       #
################################################

cmake_minimum_required( VERSION 2.5 )

## PROJECT NAME and supported LANGUAGES

project ( OpenMPTest CXX )

################################################
#   Simple C++ program                         #
################################################

find_package ( OpenMP )

message( STATUS "OPENMP_FOUND=${OPENMP_FOUND}, OpenMP_CXX_FLAGS=${OpenMP_CXX_FLAGS}" )
message( STATUS "OpenMP_C_FLAGS=${OpenMP_C_FLAGS}" )

set ( CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} )

add_executable( example example.cpp )
