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
 # @brief CMake configuration file for lama - parallel I/O examples
 # @author Thomas Brandes
 # @date 03.07.2016
###

get_filename_component ( SUBDIR ${CMAKE_CURRENT_SOURCE_DIR} NAME )

set ( EXAMPLE_FILES writeParallel.cpp readParallel.cpp utility.hpp )

add_executable( writeParallel.exe writeParallel.cpp )
add_executable( readParallel.exe readParallel.cpp )

target_link_libraries( writeParallel.exe ${PROJECT_NAME} )
target_link_libraries( readParallel.exe ${PROJECT_NAME} )

set ( INSTALL_EXAMPLE_DIR "share/examples/scai-${MODULE_NAME}-${SCAI_VERSION}/${SUBDIR}" )

include ( Functions/buildExamples )

genExampleMakefile()       #  build make.inc, makefile.install -> makefile

## install examples sources

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
