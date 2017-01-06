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

# copy run_all.sh
file( COPY run_all.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# configure make.inc for examples
include ( Functions/prepareExamplesMakeInc )
prepareExamplesMakeInc()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../examples_make.inc.in" make.inc )
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../examples_Makefile.in" install_Makefile COPYONLY )

## install examples sources, (install_)Makefile, make.inc, run_all.sh

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( PROGRAMS    run_all.sh
          DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/make.inc
          DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/install_Makefile
          DESTINATION ${PROJECT_EXAMPLE_DIR}
          RENAME      Makefile  )

