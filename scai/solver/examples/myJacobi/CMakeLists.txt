###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief CMake configuration file for scai_lama examples
 # @author Jan Ecker
 # @date 02.12.2015
###

set ( HEADER MyJacobiModule.hpp )

set ( MODULE_FILE MyJacobiModule.cpp )
set ( EXECUTABLE_FILE myJacobi.cpp )
set ( EXECUTABLE myJacobi.exe )

# build MyJacobiModule
add_library( MyJacobiModule SHARED ${MODULE_FILE} )
target_link_libraries( MyJacobiModule ${PROJECT_NAME} )

#build main for MyJacobiModule
add_executable( ${EXECUTABLE} ${EXECUTABLE_FILE} )
target_link_libraries( ${EXECUTABLE} MyJacobiModule )

set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${EXECUTABLE_FILE} ${MODULE_FILE} )
set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${EXECUTABLE}" )
set ( EXAMPLE_LIBS "MyJacobiModule" ) # for make.inc

get_filename_component ( SUBDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME )

set ( INSTALL_EXAMPLE_DIR "share/examples/scai-${MODULE_NAME}-${SCAI_VERSION}/${SUBDIR}" )

include ( Functions/buildExamples )

genExampleMakefile()       #  build make.inc, makefile.install -> makefile

## install examples sources

install ( FILES       ${EXAMPLE_FILES} ${HEADER}
          DESTINATION ${INSTALL_EXAMPLE_DIR} )
