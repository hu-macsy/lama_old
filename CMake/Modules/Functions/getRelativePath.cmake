###
 # @file getRelativePath.cmake
 #
 # @license
 # Copyright (c) 2009-2013
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief Function to return relative path between two pathes
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
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

    message ( STATUS "RELATIVE_PATH = ${${RELATIVE_PATH}} for PATH1=${PATH1} and PATH2=${PATH2}" )

endfunction ( get_relative_path )
