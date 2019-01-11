###
 # @file scai_example_directory.cmake
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
 # @brief Funciton to get the example installation directory
 # @author Thomas Brandes
 # @date 03.07.2017
###

## Function to get the prefered installation directory for examples

function ( scai_example_directory output_var )

   get_filename_component ( SUBDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME )

   if ( "${SUBDIR}" STREQUAL "examples" )
       set( ${output_var} "share/examples/${SCAI_LIBRARY_PREFIX}-${MODULE_NAME}-${SCAI_VERSION}" PARENT_SCOPE )
   else ()
       set( ${output_var} "share/examples/${SCAI_LIBRARY_PREFIX}-${MODULE_NAME}-${SCAI_VERSION}/${SUBDIR}" PARENT_SCOPE )
   endif ()

endfunction ()

