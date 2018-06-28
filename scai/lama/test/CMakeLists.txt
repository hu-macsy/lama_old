###
 # @file scai/lama/test/CMakeLists.txt
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
 # @brief scai/lama/test/CMakeLists.txt
 # @author Thomas Brandes
 # @date 25.07.2017
###

### set this dir as lama source dir ###

scai_project ( 

    CLASSES

    HEADERS
        
       TestVectors

    SOURCES

        lamaTest

        DenseVectorTest
        GridVectorTest
        GridSectionTest
        L1NormTest
        L2NormTest
        MaxNormTest
        NormTest
        ScalarTest
        PartitionIOTest
        SparseVectorTest
        VectorTest
        VectorAssemblyTest
        VectorConstructorTest
        VectorKindTest
        VersionTest

    )

### add subdirectories ###

scai_subdirectories( storage matrix )

### Some tests require path to the testfile directory

add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../res/testfiles" )

scai_test ( EXECUTABLE lamaTest 
            FILES      ${CXX_SOURCES} 
            UNIT_TEST
            RUN )
            
### Copy the README and xmltest.sh file in build directory ### 

file ( COPY README DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

###  copy xmltest.sh script

scai_test_scripts( SCRIPTS       xmltest.sh )

###  configure code_coverage.sh (if USE_CODE_COVERAGE ) in build direcotry

scai_test_code_coverage( ${USE_CODE_COVERAGE} )

