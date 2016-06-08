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

## Add lama include directory for all examples 

include_directories( "${CMAKE_CURRENT_BINARY_DIR}/../include" )

add_subdirectory ( bench )
add_subdirectory ( labelrank )
add_subdirectory ( tutorial )

file( COPY run_all.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

## install examples
install ( FILES Makefile DESTINATION ${PROJECT_EXAMPLE_DIR} )
install ( PROGRAMS run_all.sh DESTINATION ${PROJECT_EXAMPLE_DIR} )
