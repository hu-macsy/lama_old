###
 # @file scai/lama/test/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief scai/lama/test/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
 # @since 1.0.0
###

### set this dir as lama source dir ###

scai_project ( 

    CLASSES

       # Configuration
        MatrixStorageTest
        NormTest
        SparseMatrixHelper
        SparseMatrixTest
        TestContext
    
    HEADERS
        
        SameMatrixHelper
        EquationHelper
        TestMacros
        TestSparseMatrices

    SOURCES

        CommunicatorTest
        COOSparseMatrixTest
        COOStorageTest
        COOUtilsTest
        CSRSparseMatrixTest
        CSRStorageTest
        CSRUtilsTest
        DenseMatrixTest
        DenseMatrixTest1
        DenseStorageTest
        DerivedSparseMatrixTest
        DIASparseMatrixTest
        DIAStorageTest
        ELLSparseMatrixTest
        ELLStorageTest
        ELLUtilsTest
        HArrayUtilsTest
        JDSSparseMatrixTest
        JDSStorageTest
        JDSUtilsTest
        L1NormTest
        L2NormTest
        MatrixNormTest
        MatrixStorageConversionTest
        MatrixTest
        MaxNormTest
        NoCommunicatorTest
        ScalarTest
        SparseAssemblyStorageTest
        StorageIOTest
        TaskSyncTokenTest
        UtilsTest
        VectorTest
    )

### checks ####
if ( ${SCAI_ASSERT_LEVEL} STREQUAL "OFF" )
   message ( STATUS "WARNING: tests will not check for exceptions (assertions are switched off)" )
endif ( ${SCAI_ASSERT_LEVEL} STREQUAL "OFF" )

### set include directories ###
include_directories ( ${Boost_INCLUDE_DIR} )

if ( SCAI_BLAS_NAME STREQUAL "MKL" )
    include_directories ( ${MKL_INCLUDE_DIRS} )
endif ( SCAI_BLAS_NAME STREQUAL "MKL" )

### add subdirectories ###
if ( MPI_FOUND AND USE_MPI )
    add_subdirectory ( distributed )
endif ( MPI_FOUND AND USE_MPI )

if ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
 
### add definitions ###
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../res/testfiles" )

### add link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add executeable ###
add_executable ( lama_test lama_test.cpp ${CXX_SOURCES} )

### target link libraries ###
if ( WIN32 )
    target_link_libraries ( lama_test ${PROJECT_NAME} )
else ( WIN32 )
    target_link_libraries ( lama_test ${PROJECT_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} )
endif ( WIN32 )

### Copy the README and xmltest.sh file in build directory ### 
file ( COPY README DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverage.sh" @ONLY)
endif ( USE_CODE_COVERAGE )
