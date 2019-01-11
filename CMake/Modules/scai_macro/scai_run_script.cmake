###
 # @file scai_run_script.cmake
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
 # @brief Macros for copy a run script in binary directory and installation directory
 # @author Thomas Brandes
 # @date 03.07.2017
###

## This macro copies an executable file from current source dir in installatin dir

macro ( scai_run_script )

    set ( options )
    set ( oneValueArgs DESTINATION )
    set ( multiValueArgs COPY )

    cmake_parse_arguments ( scai_run_script "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( "${scai_run_script_DESTINATION}" STREQUAL "" )
        message ( FATAL_ERROR "scai_run_script: DESTINATION not specified" )
    endif ()

    ## run script might be called in build/binary directory

    file ( COPY ${scai_run_script_COPY} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

    ## but also later in installation directory

    install ( PROGRAMS    ${scai_run_script_COPY} 
              DESTINATION ${scai_run_script_DESTINATION} )

endmacro ()

