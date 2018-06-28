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
 # @brief CMake configuration file for scai_logging examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

include ( scai_macro/scai_add_example )
include ( scai_macro/scai_example_makefile )
include ( scai_macro/scai_run_script )
include ( scai_function/scai_example_directory )

# LogLevels example, here we compile also LOG_TRACE macros

set ( SCAI_LOG_LEVEL TRACE ) # also needed in prepareExamplesMakeInc for make.inc

scai_add_example ( EXECUTABLE LogLevels.exe
                   FILES      LogLevels.cpp )

set_target_properties( LogLevels.exe PROPERTIES COMPILE_FLAGS -DSCAI_LOG_LEVEL_${SCAI_LOG_LEVEL} )

# Try to find OpenMP compile and link flags for the OpenMP example
# OpenMPTest example

include ( Package/OpenMP )

if ( OPENMP_FOUND )

    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )

    scai_add_example ( EXECUTABLE LogOpenMP.exe FILES LogOpenMP.cpp )

else  ( OPENMP_FOUND )

    message ( STATUS "LogOpenMP example not built, no compiler option found" )

endif ( OPENMP_FOUND )

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# generate build file in example directory
scai_example_makefile ( DESTINATION ${INSTALL_EXAMPLE_DIR} )

# copy the script to run all examples
scai_run_script ( COPY run_all.sh DESTINATION ${INSTALL_EXAMPLE_DIR} )

## install examples sources

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${INSTALL_EXAMPLE_DIR} )

