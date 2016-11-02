###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
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

set ( EXAMPLE_FILES matrixConvert.cpp vectorConvert.cpp matrixRepartition.cpp vectorRepartition.cpp )

add_executable( matrixConvert.exe matrixConvert.cpp )
add_executable( vectorConvert.exe vectorConvert.cpp )
add_executable( matrixRepartition.exe matrixRepartition.cpp )
add_executable( vectorRepartition.exe vectorRepartition.cpp )

target_link_libraries( matrixConvert.exe ${PROJECT_NAME} )
target_link_libraries( vectorConvert.exe ${PROJECT_NAME} )
target_link_libraries( matrixRepartition.exe ${PROJECT_NAME} )
target_link_libraries( vectorRepartition.exe ${PROJECT_NAME} )

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_make.inc.in" make.inc )
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_Makefile.in" install_Makefile COPYONLY )

## install examples

set ( PROJECT_EXAMPLE_SUB_DIR ${PROJECT_EXAMPLE_DIR}/io )

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/make.inc
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/install_Makefile
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR}
          RENAME      Makefile  )
