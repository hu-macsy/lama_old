###
 # @file setCMakeQuiet.cmake
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
 # @brief CMake function for quiet cmake mode (set cmake quiet for untyped messages)
 # @author Lauretta Schubert
 # @date 01.04.2016
###

## adapted to example from: http://stackoverflow.com/questions/10509380/tell-cmake-to-be-quiet

function ( message )

    list ( GET ARGV 0 MessageType )

    if    ( MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR OR
            MessageType STREQUAL WARNING     OR MessageType STREQUAL AUTHOR_WARNING OR
            MessageType STREQUAL STATUS )

        list ( REMOVE_AT ARGV 0 )
        _message( ${MessageType} ${ARGV} )

    endif ( )

endfunction ()
