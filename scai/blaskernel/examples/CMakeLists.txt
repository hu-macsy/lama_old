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
 # @brief CMake configuration file for blaskernel examples
 # @author Eric Schricker
 # @date 28.01.2016
###

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SCAI_LANG_FLAGS}" )

if    ( OPENMP_FOUND AND USE_OPENMP )
	set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif ( OPENMP_FOUND AND USE_OPENMP )
 
# example of a library module ( will be used by executable UseModule )

foreach ( executable CG_BLAS )
    add_executable ( ${executable}.exe ${executable}.cpp )
    target_link_libraries ( ${executable}.exe ${PROJECT_NAME} )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${executable}.cpp )
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
endforeach ( executable )

# copy run_all.sh
file( COPY run_all.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

# configure make.inc for examples
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

