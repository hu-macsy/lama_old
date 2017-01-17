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
