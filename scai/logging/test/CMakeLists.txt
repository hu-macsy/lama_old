###
 # @file scai/logging/test/CMakeLists.txt
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
 # @brief CMake configuration file for tests of SCAI module logging
 # @author Thomas Brandes
 # @date 04.07.2017
###

## Note: logging has no Boost unit test, but own test script

# Build simpleLogging, one version for each possible logging level ( used by test.sh )

foreach ( level OFF FATAL ERROR WARN INFO DEBUG TRACE )
    scai_test ( EXECUTABLE simpleLogging${level} FILES simpleLogging.cpp )
    target_compile_definitions ( simpleLogging${level} PRIVATE -DSCAI_LOG_LEVEL_${level} )
endforeach ()

# Build complexLogging ( used by test.sh )

scai_test ( EXECUTABLE complexLogging FILES complexLogging.cpp )
target_compile_definitions ( complexLogging PRIVATE -DSCAI_LOG_LEVEL_TRACE )

scai_test_scripts( SCRIPTS       test.sh
                   FILES         README    )

###  configure code_coverage.sh (if USE_CODE_COVERAGE ) in build direcotry

scai_test_code_coverage( ${USE_CODE_COVERAGE} )

add_test ( logging_test test.sh WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR} )

