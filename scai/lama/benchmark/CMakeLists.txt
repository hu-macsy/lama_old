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
 # @brief CMake configuration file for benchmarks of the LAMA package
 # @author Thomas Brandes
 # @date 13.09.2017
###

## Add scai include directory for all benchmarks

include_directories( "${CMAKE_CURRENT_BINARY_DIR}/../include" )

## not yet added. ConvertMatrixStorageBenchmark.cpp

set ( CXX_SOURCES LAMAMPIBenchmark.cpp
                  PMVBenchmark.cpp
                  ConvertMatrixStorageBenchmark.cpp
                  LAMAInputSet.cpp 
                  LAMAFileInputSet.cpp 
                  LAMAInputSetComplexityVisitor.cpp
                  PoissonInputSet.cpp
                  RandomInputSet.cpp
    )

add_library ( lama_benchmark SHARED ${CXX_SOURCES} )

## Benchmarks are linked with module library and benchmark framework

target_link_libraries( lama_benchmark ${SCAI_LIBRARY_PREFIX}_benchmark ${MODULE_LIBRARY} )

## All benchmarks are installed as libraries in the benchmark directory where they might be loaded dynamically

install ( TARGETS lama_benchmark DESTINATION benchmark )

