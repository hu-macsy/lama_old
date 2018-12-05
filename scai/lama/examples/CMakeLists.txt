###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
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
add_subdirectory ( mandelbrot )
add_subdirectory ( stencil )
add_subdirectory ( nbody )
add_subdirectory ( tutorial )
add_subdirectory ( io )
add_subdirectory ( image )
add_subdirectory ( tasking )

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# copy the script to run all examples
scai_run_script ( COPY run_all.sh DESTINATION ${INSTALL_EXAMPLE_DIR} )

## copy Makefile from source directory, not configured here

install ( FILES       ${CMAKE_CURRENT_SOURCE_DIR}/Makefile 
          DESTINATION ${INSTALL_EXAMPLE_DIR} )
