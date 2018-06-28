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
 # @brief CMake configuration file for common examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

include ( scai_macro/scai_add_example )
include ( scai_macro/scai_example_makefile )
include ( scai_macro/scai_run_script )
include ( scai_function/scai_example_directory )

# example of a library module ( will be used by executable UseModule )
# must be build SHARED regardless the SCAI_LIBRARY_TYPE is SHARED or STATIC

scai_add_example ( LIBRARY Module FILES Module.cpp DynRoutine.hpp )
scai_add_example ( LIBRARY DummyModule FILES DummyModule.cpp )

set ( EXECUTABLE_FILES ${EXECUTABLE_FILES} 
          BenchMath
          BenchPointers
          DemoFactory
          DemoFactory1
          DemoGrid
          DemoLoopMacro
          DemoMath
          DemoMetaProgramming
          DemoPointer
          DemoSettings
          DemoTypeTrait
          ExceptionDemo
          TimePrecision
          UseModule
    )

if ( ${USE_COMPLEX} )
    set ( EXECUTABLE_FILES ${EXECUTABLE_FILES} DemoComplex )
endif ()

foreach ( executable ${EXECUTABLE_FILES} )

    scai_add_example( EXECUTABLE ${executable}.exe
                      FILES      ${executable}.cpp )

endforeach ()

if ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
endif ()

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# generate build file in example directory
scai_example_makefile ( DESTINATION ${INSTALL_EXAMPLE_DIR} )

# copy the script to run all examples
scai_run_script ( COPY run_all.sh DESTINATION ${INSTALL_EXAMPLE_DIR} )

install ( FILES       ${EXAMPLE_FILES} 
          DESTINATION ${INSTALL_EXAMPLE_DIR} )
