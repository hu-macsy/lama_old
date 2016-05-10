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
 # @brief CMake configuration file for taskig CUDA examples
 # @author Thomas Brandes
 # @date 08.03.2016
###

#empty EXAMPLE_FILES and EXAMPLE_EXECUTES in subdir

set ( EXAMPLE_FILES "" )
set ( EXAMPLE_EXECUTABLES "" )

string (REPLACE ";" " " CUDA_NVCC_FLAGS_CLEAN "${CUDA_NVCC_FLAGS}")

foreach ( executable BenchStream DemoTasking DemoStream DemoAsync )
    add_executable( ${executable}.exe ${executable}.cpp )
    target_link_libraries( ${executable}.exe ${PROJECT_NAME} )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${executable}.cpp ) 
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
endforeach ( executable )

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" )

## install examples
install ( FILES ${EXAMPLE_FILES} ${CMAKE_CURRENT_BINARY_DIR}/make.inc DESTINATION ${PROJECT_EXAMPLE_DIR}/cuda )

# install Makefile
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_Makefile.in"
				 "${CMAKE_INSTALL_PREFIX}/${PROJECT_EXAMPLE_DIR}/cuda/Makefile" COPYONLY )
