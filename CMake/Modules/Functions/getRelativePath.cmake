###
 # @file Functions.cmake
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
 # @brief CMake functions and macros
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

# returns the relative path to the actual directory to the CMAKE_SOURCE_DIR
function    ( get_relative_path RELATIVE_PATH )
    # get relative path
    string ( LENGTH "${CMAKE_SOURCE_DIR}" CMAKE_SOURCE_DIR_LENGTH )
    string ( LENGTH ${CMAKE_CURRENT_SOURCE_DIR} CMAKE_CURRENT_SOURCE_DIR_LENGTH )
    
    if    ( ${CMAKE_SOURCE_DIR_LENGTH} LESS ${CMAKE_CURRENT_SOURCE_DIR_LENGTH} )
        math ( EXPR CMAKE_SOURCE_DIR_LENGTH ${CMAKE_SOURCE_DIR_LENGTH}+1 )
        set ( PATH_SUFFIX / )
    endif ( ${CMAKE_SOURCE_DIR_LENGTH} LESS ${CMAKE_CURRENT_SOURCE_DIR_LENGTH} )
    
    math ( EXPR PATH_LENGTH ${CMAKE_CURRENT_SOURCE_DIR_LENGTH}-${CMAKE_SOURCE_DIR_LENGTH} )
    string ( SUBSTRING ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR_LENGTH} ${PATH_LENGTH} PATH )
    
    # set return parameter (via PARENT_SCOPE)
    set ( ${RELATIVE_PATH} ${PATH}${PATH_SUFFIX} PARENT_SCOPE )
endfunction ( get_relative_path )