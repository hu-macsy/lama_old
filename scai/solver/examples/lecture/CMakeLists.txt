###
 # @file solver/examples/lecture/CMakeLists.txt
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
 # @brief CMake configuration file for examples of LAMA lecture.
 # @author Thomas Brandes
 # @date 02.12.2017
###

include ( scai_macro/scai_add_example )
include ( scai_macro/scai_example_makefile )
include ( scai_macro/scai_run_script )

set ( EXECUTABLES task0 task1 task2 task3 task4 task4a )

if    ( CUDA_FOUND AND USE_CUDA )
    set ( EXECUTABLES ${EXECUTABLES} task5 )
endif ( CUDA_FOUND AND USE_CUDA )

foreach ( executable ${EXECUTABLES} )

    scai_add_example( EXECUTABLE ${executable}.exe
                      FILES      ${executable}.cpp )

endforeach ( executable )

# other files to copy in binary + install directory

set ( OTHER_FILES gr_30_30.mtx log_config )

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# generate build file in example directory
scai_example_makefile ( DESTINATION ${INSTALL_EXAMPLE_DIR} )

file ( COPY ${OTHER_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

## install examples sources + data file
install ( FILES       ${EXAMPLE_FILES} ${OTHER_FILES}
          DESTINATION ${INSTALL_EXAMPLE_DIR} )
