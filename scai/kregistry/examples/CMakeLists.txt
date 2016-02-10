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
 # @brief CMake configuration file for scai_kregistry examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

foreach ( executable DemoKernelInterface BenchInterface )
    add_executable( ${executable} ${executable}.cpp )
    target_link_libraries( ${executable} ${PROJECT_NAME} )
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
endforeach ( executable )

## install examples
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )
install ( DIRECTORY ./ DESTINATION ${PROJECT_EXAMPLE_DIR} FILES_MATCHING REGEX "^.*.(c|h)pp$" )
install ( FILES Makefile ${CMAKE_CURRENT_BINARY_DIR}/make.inc DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( PROGRAMS run_all.sh DESTINATION ${PROJECT_EXAMPLE_DIR} )