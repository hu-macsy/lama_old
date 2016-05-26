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

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_make.inc.in" make.inc )
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_Makefile.in" install_Makefile COPYONLY )

## install examples

set ( PROJECT_EXAMPLE_SUB_DIR "${PROJECT_EXAMPLE_DIR}/myJacobi" )

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/make.inc
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/install_Makefile
          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR}
          RENAME      Makefile  )

