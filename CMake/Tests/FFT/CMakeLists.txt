###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
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
 # @brief CMakeLists.txt
 # @author Thomas Brandes
 # @date 25.04.2013
###

################################################
#   CMake example for C++ program using FFTW   #
################################################

cmake_minimum_required( VERSION 2.8 )

## PROJECT NAME and supported LANGUAGES

project ( FFTTest CXX )

set ( CMAKE_MODULE_PATH "../../Modules" )

find_package( FFTW )

if ( NOT FFTW_FOUND )
    message ( STATUS "CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID}" )
    find_package ( MKL )
    if ( MKL_FOUND )
        find_path( FFTW_INCLUDE_DIR fftw3.h
                   ${MKL_INCLUDE_DIR}/fftw )
        if ( FFTW_INCLUDE_DIR )
            set ( FFTW_FOUND TRUE )
            set ( FFTW_LIBRARIES ${MKL_LIBRARIES} )
        endif ()
    endif ()
endif ()

message ( STATUS "FFTW_FOUND=${FFTW_FOUND}" )
message ( STATUS "FFTW_INCLUDE_DIR=${FFTW_INCLUDE_DIR}" )
message ( STATUS "FFTW_LIBRARIES=${FFTW_LIBRARIES}" )

include_directories( ${FFTW_INCLUDE_DIR} )
add_executable( fft_simple fft_simple.cpp )
target_link_libraries( fft_simple ${FFTW_LIBRARIES} )

