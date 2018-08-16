###
 # @file scai_macro/scai_build_variable.cmake
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
 # @brief Macros for building examples in SCAI module projects
 # @author Thomas Brandes
 # @date 03.07.2017
###

include ( scai_function/listToString )
include ( CMakeParseArguments )

## Macro to define SCAI variables ( always in cache )

macro ( scai_build_variable )

    set ( options BOOL )
    set ( oneValueArgs NAME DEFAULT DOCSTRING )
    set ( multiValueArgs CHOICES )

    cmake_parse_arguments ( scai_build_variable "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( DEFINED ${scai_build_variable_NAME} )
        set ( _value ${${scai_build_variable_NAME}} )
    else  ()
        set ( _value "" )
    endif ()

    # message ( STATUS "var ${scai_build_variable_NAME} = \"${_value}\"" )

    if ( "${_value}" STREQUAL "" )
        set ( _value ${scai_build_variable_DEFAULT} )
    endif ()

    if ( ${scai_build_variable_BOOL} )

        ## make sure that _value is either ON or OFF

        if ( _value )
            set ( _value ON )
        else  ()
            set ( _value OFF )
        endif ()

       #  FORCE it in cache so variables can be redefined via cmake --Dscai_var=value

        set( ${scai_build_variable_NAME} ${_value} CACHE BOOL "Enable / Disable ${scai_build_variable_DOCSTRING}" FORCE )

    else ()

        ## make sure that _value is a legel value

        listToString ( ", " "${scai_build_variable_CHOICES}" _str_choices )

        # message ( STATUS "check legal value: \"${_value}\" must be in ${_str_choices}" )

        list ( FIND scai_build_variable_CHOICES "${_value}" _index )

        if ( ${_index} EQUAL -1 )
            message ( FATAL_ERROR "${scai_build_variable_NAME}: value ${_value} illegal, choices are ${_str_choices}" )
        endif ()

       #  FORCE it in cache so variables can be redefined via cmake --Dscai_var=value

        set( ${scai_build_variable_NAME} ${_value} CACHE STRING
             "${scai_build_variable_DOCSTRING} ${_str_choices}" FORCE )

    endif ()

endmacro ()
