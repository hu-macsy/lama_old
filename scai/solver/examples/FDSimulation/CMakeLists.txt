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
 # @author Vanessa Wolff
 # @date 24.05.2016
###

set ( EXECUTABLE_FILE FDSimulation.cpp )
set ( EXECUTABLE FDSimulation.exe )

#build main for FDSimulation
add_executable( ${EXECUTABLE} ${EXECUTABLE_FILE} )
target_link_libraries( ${EXECUTABLE} ${PROJECT_NAME})

set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${EXECUTABLE_FILE} ${MODULE_FILE} )
set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${EXECUTABLE}" )

# seismogram.gnuplot: copy only to build directory
configure_file( seismogram.gnuplot seismogram.gnuplot COPYONLY )
# Configuration.txt: copy only to build directory
configure_file( Configuration.txt Configuration.txt COPYONLY )

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )

#install examples
install ( FILES ${EXAMPLE_FILES} ${HEADER} ${CMAKE_CURRENT_BINARY_DIR}/make.inc seismogram.gnuplot Configuration.txt Configuration.hpp
		  DESTINATION ${PROJECT_EXAMPLE_DIR}/FDSimulation )

# install Makefile
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_Makefile.in"
				 "${CMAKE_INSTALL_PREFIX}/${PROJECT_EXAMPLE_DIR}/FDSimulation/Makefile" COPYONLY )