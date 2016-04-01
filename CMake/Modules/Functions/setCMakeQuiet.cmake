###
 # @file setCMakeQuiet.cmake
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
 # @brief CMake function for quiet cmake mode (set cmake quiet for untyped messages)
 # @author Lauretta Schubert
 # @date 01.04.2016
 # @since 2.0.0
###

## adapted to example from: http://stackoverflow.com/questions/10509380/tell-cmake-to-be-quiet

function    ( message )
	list ( GET ARGV 0 MessageType )
  	if    ( MessageType STREQUAL FATAL_ERROR OR MessageType STREQUAL SEND_ERROR OR
      	    MessageType STREQUAL WARNING     OR MessageType STREQUAL AUTHOR_WARNING OR
	      	MessageType STREQUAL STATUS )
  		list ( REMOVE_AT ARGV 0 )
    	_message( ${MessageType} ${ARGV} )
  	endif ( )
endfunction ( message )