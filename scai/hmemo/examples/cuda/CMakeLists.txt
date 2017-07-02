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
 # @brief CMake configuration file for common examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

# REMARKS
#
#   - include_directories are taken over from project
#

#empty EXAMPLE_FILES and EXAMPLE_EXECUTES in subdir

set ( EXAMPLE_FILES "" )
set ( EXAMPLE_EXECUTABLES "" )

string (REPLACE ";" " " CUDA_NVCC_FLAGS_CLEAN "${CUDA_NVCC_FLAGS}")

foreach ( executable CUSparseExample CUBlasExample MemBandwidth Prefetch Allocate )
    add_executable( ${executable}.exe ${executable}.cpp )
    target_link_libraries( ${executable}.exe ${PROJECT_NAME} )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${executable}.cpp ) 
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
endforeach ( executable )

foreach ( executable Example1 Example2 Devices CUDABenchContext AliasProblem )
    cuda_add_executable( ${executable}.exe ${executable}.cu )
    target_link_libraries( ${executable}.exe ${PROJECT_NAME} )
    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${executable}.cu )
    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
endforeach ( executable )

set ( INSTALL_EXAMPLE_DIR "share/examples/scai-${MODULE_NAME}-${SCAI_VERSION}/cuda" )

include ( Functions/buildExamples )

genExampleMakefile()       #  build make.inc, makefile.install -> makefile

## install examples sources

install ( FILES       ${EXAMPLE_FILES}
          DESTINATION ${INSTALL_EXAMPLE_DIR} )
