###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake configuration file for common examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

# ThreadTest example

include ( Compiler/CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 )
    include( Package/Boost )
    include_directories ( ${BOOST_INCLUDE_DIR} )
endif ( NOT CXX_SUPPORTS_C11 )

if ( ${SCAI_LIBRARY_TYPE} MATCHES SHARED )
    # example of a library module ( will be used by executable UseModule )
    add_library( Module SHARED Module.cpp )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} Module.cpp )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} DynRoutine.hpp )
    target_link_libraries( Module ${PROJECT_NAME} )
endif ( ${SCAI_LIBRARY_TYPE} MATCHES SHARED )

# example of a library module ( no use of other LAMA stuff )
add_library( DummyModule SHARED DummyModule.cpp )
set ( EXAMPLE_FILES ${EXAMPLE_FILES} DummyModule.cpp )

set ( EXAMPLE_MODULES "Module.so DummyModule.so" ) # for make.inc

foreach ( executable ExceptionDemo CriticalRegion Barrier 
          DemoFactory DemoFactory1 DemoPointer UseModule TimePrecision 
          DemoComplex DemoTypeTrait DemoSettings DemoMath DemoFunction 
          DemoMetaProgramming BenchPointers BenchMath )
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

# install Makefile
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../examples_Makefile.in"
				 "${CMAKE_INSTALL_PREFIX}/${PROJECT_EXAMPLE_DIR}/Makefile" COPYONLY )
