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
 # @brief CMake configuration file for scai_logging examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

# Try to find OpenMP compile and link flags for the OpenMP example

# LogLevels example, here we compile also LOG_TRACE macros
set ( SCAI_LOGGING_LEVEL TRACE ) # also needed in prepareExamplesMakeInc for make.inc

add_executable( LogLevels.exe LogLevels.cpp )
set ( EXAMPLE_FILES ${EXAMPLE_FILES} LogLevels.cpp )
set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} LogLevels.exe" )
target_link_libraries( LogLevels.exe ${PROJECT_NAME} )
set_target_properties( LogLevels.exe PROPERTIES COMPILE_FLAGS -DSCAI_LOG_LEVEL_${SCAI_LOGGING_LEVEL} )

# ThreadTest example, here we need the Boost Thread library

# Boost thread required

set ( Boost_DEBUG TRUE )
set ( Boost_DETAILED_FAILURE_MSG TRUE )

find_package ( Boost ${SCAI_FIND_PACKAGE_FLAGS} COMPONENTS thread system )

# OpenMPTest example

include ( Package/OpenMP )

if    ( OPENMP_FOUND )
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # Adding this example only possible if OpenMP is available
    add_executable( LogOpenMP.exe LogOpenMP.cpp )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} LogOpenMP.cpp )
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} LogOpenMP.exe" )
    target_link_libraries( LogOpenMP.exe ${PROJECT_NAME} )
else  ( OPENMP_FOUND )
    message ( STATUS "LogOpenMP example not built, no compiler option found" )
endif ( OPENMP_FOUND )

file( COPY run_all.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../examples_make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )

## install examples
install ( FILES ${EXAMPLE_FILES} DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( FILES config Makefile ${CMAKE_CURRENT_BINARY_DIR}/make.inc DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( PROGRAMS run_all.sh DESTINATION ${PROJECT_EXAMPLE_DIR} )
