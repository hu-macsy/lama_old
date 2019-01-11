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

include ( scai_macro/scai_add_example )
include ( scai_macro/scai_example_makefile )
include ( scai_macro/scai_run_script )

foreach ( executable labelrank )

    scai_add_example( EXECUTABLE ${executable}.exe
                      FILES      ${executable}.cpp )

endforeach ( executable )

# add config.hpp for installation
set ( EXAMPLE_FILES ${EXAMPLE_FILES} Config.hpp )

# add mtx files for installation
set ( EXAMPLE_FILES ${EXAMPLE_FILES} affinity.mtx labels.mtx )

# copy mtx files to build dir
file( COPY affinity.mtx labels.mtx DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# generate build file in example directory
scai_example_makefile ( DESTINATION ${INSTALL_EXAMPLE_DIR} )

## install examples sources

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${INSTALL_EXAMPLE_DIR} )

