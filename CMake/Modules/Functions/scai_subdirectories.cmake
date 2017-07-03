###
 # @file scai_subdirectories.cmake
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
 # @brief CMake macro to define a SCAI module library                          
 # @author Thomas Brandes
 # @date 04.07.2017
###

## This macro adds subdirectories specified in the argument list with
## some special handling:
##
##   test      - only if TEST enabled
##   cuda      - only if CUDA enabled
##   examples  - only if BUILD_EXAMPLES has added it
##   mic       - only if MIC is used

macro ( scai_subdirectories )

    foreach ( dir ${ARGN} )

        ## redundant to check if dir exists as CMake gives already message

        if ( ${dir} STREQUAL "examples" )

            if ( BUILD_EXAMPLES )
                add_subdirectory ( examples )
            endif ()

        elseif ( ${dir} STREQUAL "test" )

            if ( FOUND_BOOST_TEST AND BUILD_TEST )
                add_subdirectory ( test )
            endif ()

        elseif ( ${dir} STREQUAL "cuda" )

            if ( CUDA_FOUND AND USE_CUDA )
                add_subdirectory ( cuda )
                cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
                set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
            endif ()

        elseif ( ${dir} STREQUAL "mic" )

            if ( USE_MIC )
                add_subdirectory ( mic )
            endif ()

        else ()

            add_subdirectory( ${dir} )

        endif ()

    endforeach ()

endmacro ()
