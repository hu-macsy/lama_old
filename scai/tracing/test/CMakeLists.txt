###
 # @file scai/tracing/test/CMakeLists.txt
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
 # @brief scai/tracing/CMakeLists.txt
 # @author Thomas Brandes
 # @date 10.06.2015
###


# Build simpleTracingXXX, one version with Tracing enabled, one without

scai_test ( EXECUTABLE simpleTracingON FILES simpleTracing.cpp )
target_compile_definitions ( simpleTracingON PRIVATE -DSCAI_TRACE_ON )

scai_test ( EXECUTABLE simpleTracingON0 FILES simpleTracing.cpp )
target_compile_definitions ( simpleTracingON0 PRIVATE -DSCAI_TRACE_ON -DUNNAMED_THREADS )

scai_test ( EXECUTABLE simpleTracingOFF FILES simpleTracing.cpp )
target_compile_definitions ( simpleTracingOFF PRIVATE -DSCAI_TRACE_OFF )

# The shell script test.sh uses NUM_THREADS_TEST to verify tracing for threads

if ( USE_OPENMP )
   set ( NUM_THREADS_TEST 4 )
else ()
   set ( NUM_THREADS_TEST 1 )
endif ()

## make.inc and test.sh need configuration

scai_test_scripts ( SCRIPTS         test.sh
                    CONFIGURE 
)

## other files are just copied

scai_test_scripts ( FILES           README )

###  configure code_coverage.sh (if USE_CODE_COVERAGE ) in build direcotry

scai_test_code_coverage( ${USE_CODE_COVERAGE} )

add_test ( tracing_test test.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )
