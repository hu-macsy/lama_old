###
 # @file scai_summary.cmake
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

include ( CMakeParseArguments )
include ( scai_function/scaiMessages )
include ( scai_function/listToString )

## Macro adds a summary string (one line) to global variable SCAI_SUMMARY

macro ( scai_summary )

    set( SCAI_SUMMARY "${SCAI_SUMMARY}\n-- ${ARGN}" )

    get_directory_property ( hasParent PARENT_DIRECTORY )

    if ( hasParent )
        set ( SCAI_SUMMARY ${SCAI_SUMMARY} PARENT_SCOPE )
    endif ()

endmacro ()

macro ( scai_summary_enabled )

    set ( options )
    set ( oneValueArgs NAME ENABLED )
    set ( multiValueArgs )

    cmake_parse_arguments ( scai_summary_enabled "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( ${scai_summary_enabled_ENABLED} )

        set ( EnabledString "${TextGreen}ENABLED${TextColorReset}" )

    else ()

        set ( EnabledString "${TextAmber}DISABLED${TextColorReset}" )

    endif ()

    scai_summary ( "     ${scai_summary_enabled_NAME} ${EnabledString}" )

endmacro ()

macro ( scai_summary_newline )

    scai_summary ( " " )

endmacro ()

#  scai_summary_external ( NAME      MPI 
#                          FOUND     ${MPI_FOUND} 
#                          VERSION   ${MPI_VERSION} 
#                          INCLUDE   ${SCAI_MPI_INCLUDE_DIR} 
#                          LIBRARIES ${SCAI_MPI_LIBRARIES} )
#
#  scai_summary_external ( NAME      OpenMP
#                          FOUND     ${OPENMP_FOUND} 
#                          VERSION   ${OPENMP_VERSION} 
#                          CXX_FLAGS ${OpenMPCXXFlags} )

macro ( scai_summary_external )

    set ( options )
    set ( oneValueArgs NAME FOUND VERSION CXX_FLAGS ENABLED )
    set ( multiValueArgs INCLUDE LIBRARIES EXECUTABLE )

    cmake_parse_arguments ( scai_summary_external "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( DEFINED scai_summary_external_ENABLED )

        scai_summary_enabled( NAME ${scai_summary_external_NAME} ENABLED ${scai_summary_external_ENABLED} )

    endif ()

    if ( ${scai_summary_external_FOUND} )

         set ( FoundString "${TextGreen}FOUND${TextColorReset}" )

         # if ( ${scai_summary_external_VERSION} )
         if ( DEFINED scai_summary_external_VERSION )
            set ( FoundString "${FoundString}      : Version ${scai_summary_external_VERSION}" )
         endif ()

    else ()

        set ( FoundString "${TextAmber}NOT FOUND${TextColorReset}" )

    endif ()

    scai_summary ( "         ${FoundString}" )

    if ( ${scai_summary_external_FOUND} )

        if ( DEFINED scai_summary_external_CXX_FLAGS )
            scai_summary ( "         CXX Flags  : ${scai_summary_external_CXX_FLAGS}" )
        endif ()

        list ( LENGTH scai_summary_external_EXECUTABLE N_EXECUTABLE )

        if ( N_EXECUTABLE GREATER 0 )
            listToString ( ", " "${scai_summary_external_EXECUTABLE}" EXECUTABLE_LIST )
            scai_summary ( "         Executable : ${EXECUTABLE_LIST}" )
        endif ()

        list ( LENGTH scai_summary_external_INCLUDE N_INCLUDE )

        if ( N_INCLUDE GREATER 0 )
            listToString ( ", " "${scai_summary_external_INCLUDE}" INCLUDE_LIST )
            scai_summary ( "         Include(${N_INCLUDE}) : ${INCLUDE_LIST}" )
        endif ()

        list ( LENGTH scai_summary_external_LIBRARIES N_LIBS )

        if ( N_LIBS GREATER 0 )
            listToString ( ", " "${scai_summary_external_LIBRARIES}" LIB_LIST )
            scai_summary ( "         Libs(${N_LIBS})    : ${LIB_LIST}" )
        endif ()

    endif ()

    ## make a final new line

    scai_summary ( " " )

endmacro ()


