###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @brief CMake configuration file for scai_lama examples
 # @author Jan Ecker
 # @date 02.12.2015
###

## Add lama include directory for all examples 

include_directories( "${CMAKE_CURRENT_BINARY_DIR}/../include" )

add_subdirectory ( bench )
add_subdirectory ( labelrank )
add_subdirectory ( lecture )
add_subdirectory ( solver )
add_subdirectory ( spy )
add_subdirectory ( tutorial )

## install examples
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )
install ( DIRECTORY ./ DESTINATION share/examples/scai-lama-${SCAI_LAMA_VERSION} FILES_MATCHING REGEX "^(.*.(cpp|hpp|mtx))|(Makefile)$" )
install ( FILES Makefile ${CMAKE_CURRENT_BINARY_DIR}/make.inc DESTINATION share/examples/scai-lama-${SCAI_LAMA_VERSION} )
