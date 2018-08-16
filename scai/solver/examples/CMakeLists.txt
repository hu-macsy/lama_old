###
 # @file solver/examples/CMakeLists.txt
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
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 25.01.2016
###

include ( scai_macro/scai_add_example )
include ( scai_macro/scai_example_makefile )
include ( scai_macro/scai_run_script )
include ( scai_function/scai_example_directory )

add_subdirectory ( FDSimulation )
add_subdirectory ( lecture )
if ( ${SCAI_LIBRARY_TYPE} MATCHES SHARED )
    add_subdirectory ( myJacobi )
endif ( ${SCAI_LIBRARY_TYPE} MATCHES SHARED )
add_subdirectory ( pardiso )
add_subdirectory ( solver )
add_subdirectory ( spy )
add_subdirectory ( eigenvalue )
add_subdirectory ( leastsquare )

# use the recommended installation directory for this example directory
scai_example_directory ( INSTALL_EXAMPLE_DIR )

# copy the script to run all examples
scai_run_script ( COPY run_all.sh DESTINATION ${INSTALL_EXAMPLE_DIR} )

## copy the Makefile from the source directory to installation

install ( FILES Makefile DESTINATION ${INSTALL_EXAMPLE_DIR} )
