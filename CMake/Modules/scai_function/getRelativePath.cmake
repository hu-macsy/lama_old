###
 # @file getRelativePath.cmake
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
 # @brief Function to return relative path between two pathes
 # @author Jan Ecker
 # @date 25.04.2013
###

# returns the relative path of PATH2 to PATH1, PATH1 must be substring of PATH2 

# get_relative_path ( RELATIVE_PATH /home /home/ecker ) -> ecker
# get_relative_path ( RELATIVE_PATH /home/scai /home/scai/user/ecker ) -> user/ecker

function ( get_relative_path RELATIVE_PATH PATH1 PATH2 )

    # get relative path between PATH1 and PATH2, PATH1 must be substring von PATH2

    string ( LENGTH ${PATH1} PATH1_LENGTH )
    string ( LENGTH ${PATH2} PATH2_LENGTH )
   
    # make sure that PATH1 is substring of PATH2

    if ( PATH1_LENGTH GREATER PATH2_LENGTH )
       message ( FATAL_ERROR "get_relative_path, PATH1 = ${PATH1} not substring of PATH2 = ${PATH2}" )
    else ()
       string ( SUBSTRING ${PATH2} 0 ${PATH1_LENGTH} CHECKPATH )

       if ( CHECKPATH STREQUAL PATH1 )
          # that is okay
       else ()
          message ( FATAL_ERROR "get_relative_path, PATH1 = ${PATH1} not substring of PATH2 = ${PATH2}" )
       endif ()

    endif ( PATH1_LENGTH GREATER PATH2_LENGTH )

    if ( ${PATH1_LENGTH} LESS ${PATH2_LENGTH} )
        math ( EXPR PATH1_LENGTH ${PATH1_LENGTH}+1 )
        set ( PATH_SUFFIX / )
    endif ( )

    math ( EXPR REL_PATH_LENGTH ${PATH2_LENGTH}-${PATH1_LENGTH} )
    string ( SUBSTRING ${PATH2} ${PATH1_LENGTH} ${REL_PATH_LENGTH} PATH )
    
    # set return parameter (via PARENT_SCOPE)
    set ( ${RELATIVE_PATH} ${PATH}${PATH_SUFFIX} PARENT_SCOPE )

    # message ( STATUS "RELATIVE_PATH = ${${RELATIVE_PATH}} for PATH1=${PATH1} and PATH2=${PATH2}" )

endfunction ( get_relative_path )
