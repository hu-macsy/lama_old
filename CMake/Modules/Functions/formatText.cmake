###
 # @file formatText.cmake
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
 # @endlicense
 #
 # @brief CMake function to format text
 # @author Jan Ecker
 # @date 25.04.2013
###

# prints colored text messages
# inspired by soci colormsg function
function    ( formatText )
    # first arg is name of "return variable"
    list ( GET ARGV 0 RESULT_NAME )
    list ( REMOVE_AT ARGV 0 )

    include ( Settings/bashFormats )

    set ( coloron FALSE )
    set ( str "" )
    foreach    ( arg ${ARGV} )

        if    ( DEFINED ${arg} )
            if    ( CMAKE_COLOR_MAKEFILE )
                set ( str "${str}${${arg}}" )
                set ( coloron TRUE )
            endif ( CMAKE_COLOR_MAKEFILE )
        else  ( DEFINED ${arg} )
            set ( str "${str}${arg}" )
            if    ( coloron )
                set ( str "${str}${TextColorReset}" )
                set ( coloron FALSE )
            endif ( coloron )
            set ( str "${str} " )
        endif ( DEFINED ${arg} )

    endforeach ( arg ${ARGV} )

    if    ( CMAKE_COLOR_MAKEFILE )
        set ( str "${str}${TextReset}${BGReset}" )
    endif ( CMAKE_COLOR_MAKEFILE )  
    set ( ${RESULT_NAME} ${str} PARENT_SCOPE )
endfunction ( formatText )
