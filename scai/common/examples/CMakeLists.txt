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
 # @brief CMake configuration file for common examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

# ThreadTest example

include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 )
    include( Package/Boost )
    include_directories ( ${BOOST_INCLUDE_DIR} )
endif ( NOT CXX_SUPPORTS_C11 )

# example of a library module ( will be used by executable UseModule )
add_library( Module SHARED Module.cpp )
set ( EXAMPLE_FILES ${EXAMPLE_FILES} Module.cpp )
set ( EXAMPLE_FILES ${EXAMPLE_FILES} DynRoutine.hpp )
target_link_libraries( Module ${PROJECT_NAME} )

# example of a library module ( no use of other LAMA stuff )
add_library( DummyModule SHARED DummyModule.cpp )
set ( EXAMPLE_FILES ${EXAMPLE_FILES} DummyModule.cpp )

set ( EXAMPLE_MODULES "Module.so DummyModule.so" ) # for make.inc

foreach ( executable ExceptionDemo CriticalRegion Barrier 
          DemoFactory DemoFactory1 DemoPointer UseModule TimePrecision 
          DemoComplex DemoTypeTrait DemoSettings DemoMath DemoFunction 
          DemoMetaProgramming BenchPointers )
    add_executable ( ${executable}.exe ${executable}.cpp )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${executable}.cpp )
    target_link_libraries ( ${executable}.exe ${PROJECT_NAME} )
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
endforeach ( executable )

if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA )

# copy run_all.sh
file( COPY run_all.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../examples_make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )

## install examples
install ( FILES ${EXAMPLE_FILES} ${CMAKE_CURRENT_BINARY_DIR}/make.inc DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( PROGRAMS run_all.sh DESTINATION ${PROJECT_EXAMPLE_DIR} )

message ( STATUS "CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}" )
message ( STATUS "PROJECT_EXAMPLE_DIR ${PROJECT_EXAMPLE_DIR}" )

# install Makefile
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../examples_Makefile.in"
				 "${CMAKE_INSTALL_PREFIX}/${PROJECT_EXAMPLE_DIR}/Makefile" COPYONLY )
