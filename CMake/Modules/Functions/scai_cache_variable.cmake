###
 # @file scai_cache_variable.cmake  
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
 # @brief Macros for building examples in SCAI module projects
 # @author Thomas Brandes
 # @date 03.07.2017
###

## include ( Functions/listToString )

## Macro to define SCAI variables ( always in cache )

macro ( scai_cache_variable )

    set ( options BOOL )
    set ( oneValueArgs NAME DEFAULT DOCSTRING )
    set ( multiValueArgs CHOICES )

    cmake_parse_arguments ( scai_cache_variable "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( DEFINED ${scai_cache_variable_NAME} )
        set ( _value ${${scai_cache_variable_NAME}} )
    else  ()
        set ( _value ${scai_cache_variable_DEFAULT} )
    endif ()

    if ( ${scai_cache_variable_BOOL} )

        ## make sure that _value is either ON or OFF

        if ( _value )
            set ( _value ON )
        else  ()
            set ( _value OFF )
        endif ()

        set( ${scai_cache_variable_NAME} ${_value} CACHE BOOL ${scai_cache_variable_DOCSTRING} FORCE )

    else ()

        ## make sure that _value is a legel value

        list ( FIND scai_cache_variable_CHOICES ${_value} _index )

        if ( ${_index} EQUAL -1 )
            message ( FATAL_ERROR "Value ${_value} is not choice out of ${scai_cache_variable_CHOICES}" )
        endif ()

        listToString ( ", " scai_cache_variable_CHOICES _str_choices )

        set( ${scai_cache_variable_NAME} ${_value} CACHE STRING
             "${scai_cache_variable_DOCSTRING} ${_str_choices}" FORCE )

    endif ()

endmacro ()

