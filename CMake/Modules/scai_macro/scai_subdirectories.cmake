###
 # @file scai_subdirectories.cmake
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
 # @brief CMake macro to define a SCAI module library                          
 # @author Thomas Brandes
 # @date 04.07.2017
###

## This macro adds subdirectories specified in the argument list with
## some special handling:
##
##   TEST      - add test only if TEST enabled and BOOST_TEST available
##   CUDA      - add cuda only if CUDA enabled

macro ( scai_subdirectories )

    ## Note: no need to check if dir exists as CMake will give error message

    set ( options TEST CUDA EXAMPLES )
    set ( oneValueArgs )
    set ( multiValueArgs )

    cmake_parse_arguments ( scai_subdirectories "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # message ( STATUS "subdirs = ${scai_subdirectories_UNPARSED_ARGUMENTS}" )

    if ( ${scai_subdirectories_EXAMPLES} )
        add_subdirectory ( examples )
    endif ()

    if ( ${scai_subdirectories_TEST} )
        if ( FOUND_BOOST_TEST AND USE_BOOST_TEST )
            add_subdirectory ( test )
        endif ()
    endif ()

    if ( ${scai_subdirectories_CUDA} )
        if ( CUDA_FOUND AND USE_CUDA )
            add_subdirectory ( cuda )
        endif ()
    endif ()

    # all other subdirectories are just added in the usual way

    foreach ( dir ${scai_subdirectories_UNPARSED_ARGUMENTS} )
        add_subdirectory( ${dir} )
    endforeach ()

endmacro ()
