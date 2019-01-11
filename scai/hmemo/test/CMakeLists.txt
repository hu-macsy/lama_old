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
 # @brief CMake configuration file for SCAI module hmemo
 # @author Thomas Brandes
 # @date 04.07.2017
###

set ( CXX_SOURCES
      hmemoTest
      ContextTest  
      AccessTest  
      HArrayTest  
      HArrayRefTest  
      MemoryTest
) 

scai_test ( EXECUTABLE hmemoTest 
            FILES      ${CXX_SOURCES} 
            UNIT_TEST
            RUN )

scai_subdirectories( CUDA )

### Copy the README and xmltest.sh file in build directory ### 

file ( COPY README DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

###  copy xmltest.sh script

scai_test_scripts( SCRIPTS       xmltest.sh )

###  configure code_coverage.sh (if USE_CODE_COVERAGE ) in build direcotry

scai_test_code_coverage( ${USE_CODE_COVERAGE} )

