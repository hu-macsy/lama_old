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

project( CUDADevices )

# make sure that the default is a RELEASE

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

# Find CUDA and set flags

include( FindCUDA )

cuda_add_executable( CUDADevices main.cu )

