###
 # @file scai/lama/test/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
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
 # @brief scai/lama/test/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
###

### set this dir as lama source dir ###

scai_project ( 

    CLASSES

       # Configuration
        SparseMatrixHelper
        SparseMatrixTest
    
    HEADERS
        
        SameMatrixHelper
        EquationHelper
        TestMacros
        TestSparseMatrices

    SOURCES

        COOSparseMatrixTest
        CSRSparseMatrixTest
        DenseMatrixTest
        DenseMatrixTest1
        DerivedSparseMatrixTest
        DIASparseMatrixTest
        ELLSparseMatrixTest
        IOStreamTest
        JDSSparseMatrixTest
        L1NormTest
        L2NormTest
        MatrixTest
        MaxNormTest
        NormTest
        ScalarTest
        FileIOTest
        PartitionIOTest
        VectorTest
        VersionTest
    )

### checks ####
if ( ${SCAI_ASSERT_LEVEL} STREQUAL "OFF" )
   message ( STATUS "WARNING: tests will not check for exceptions (assertions are switched off)" )
endif ( ${SCAI_ASSERT_LEVEL} STREQUAL "OFF" )

### set include directories ###
include_directories ( ${SCAI_BOOST_INCLUDE_DIR} )

if ( SCAI_BLAS_NAME STREQUAL "MKL" )
    include_directories ( ${MKL_INCLUDE_DIRS} )
endif ( SCAI_BLAS_NAME STREQUAL "MKL" )

### add subdirectories ###

add_subdirectory ( storage )
add_subdirectory ( matrix )

### add definitions ###
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../res/testfiles" )

### add link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add executeable ###
add_executable ( lamaTest lamaTest.cpp ${CXX_SOURCES} )

### target link libraries ###
if ( WIN32 )
    target_link_libraries ( lamaTest ${PROJECT_NAME} )
else ( WIN32 )
    target_link_libraries ( lamaTest ${PROJECT_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
endif ( WIN32 )

### Copy the README and xmltest.sh file in build directory ### 
file ( COPY README DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/code_coverage.sh" @ONLY)
    file ( COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/code_coverage.sh 
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
    file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../scai_code_coverage_functions.sh 
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
endif ( USE_CODE_COVERAGE )
