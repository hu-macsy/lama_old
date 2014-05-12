###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief test/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
 # @since 1.0.0
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### Add classes and header files ###
set ( CLASSES
        Configuration
        Context
        MatrixStorageTest
        NormTest
        SameMatrixHelper
        SparseMatrixHelper
        SparseMatrixTest
    )
    
set ( HEADERS
        EquationHelper
        TestMacros
        TestSparseMatrices
    )

set ( SOURCE_FILES
        #BiCGstabTest
        BiCGTest
        BLAS1Test
		BLAS2Test
		BLAS3Test
        CGTest
        CommonLoggerTest
        CommunicatorTest
        COOSparseMatrixTest
        COOStorageTest
        COOUtilsTest
        CriterionTest
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
        ExceptionTest
        FileLoggerTest
        GMRESTest
        HostReadAccessTest
        HostWriteAccessTest
        InverseSolverTest
        IterationCountTest
        JacobiTest
        JDSSparseMatrixTest
        JDSStorageTest
        JDSUtilsTest
        L1NormTest
        L2NormTest
        LAMAArrayTest
        LAPACKTest
        MatrixStorageConversionTest
        MatrixTest
        MaxNormTest
        NoCommunicatorTest
        ResidualStagnationTest
        ResidualThresholdTest
        ScalarTest
        SimpleAMGTest
        SolutionProxyTest
        SORTest
        SparseAssemblyStorageTest
        SpecializedJacobiTest
        StorageIOTest
        TaskSyncTokenTest
        ThreadPoolTest
        TimerTest
        TracingTest
        UtilsTest
        VectorTest
    )

lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )
lama_sources ( ${SOURCE_FILES} )

### checks ####
if ( ${LAMA_ASSERT_LEVEL} STREQUAL "OFF" )
   message ( STATUS "WARNING: tests will not check for exceptions (assertions are switched off)" )
endif ( ${LAMA_ASSERT_LEVEL} STREQUAL "OFF" )

if ( NOT ( Boost_UNIT_TEST_FRAMEWORK_LIBRARY OR Boost_REGEX ) )
    message ( FATAL_ERROR "No Boost unit test available." )
endif ( NOT ( Boost_UNIT_TEST_FRAMEWORK_LIBRARY OR Boost_REGEX ) )

### set include directories ###
include_directories ( .. ../lama ${Boost_INCLUDE_DIRS} )

### add subdirectories ###
if ( MPI_FOUND AND LAMA_USE_MPI )
    add_subdirectory ( distributed )
endif ( MPI_FOUND AND LAMA_USE_MPI )

if ( CUDA_FOUND AND LAMA_USE_CUDA AND CUDA_HAVE_GPU )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND LAMA_USE_CUDA AND CUDA_HAVE_GPU )
 
### add definitions ###
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_SOURCE_DIR}/../res/testfiles" )

### add link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add executeable ###
add_executable ( lama_test lama_test.cpp ${CXX_SOURCES} )

### target link libraries ###
if ( WIN32 )
    target_link_libraries ( lama_test ama )
else ( WIN32 )
    if ( LAMA_LIBRARY_TYPE STREQUAL "STATIC" )
        target_link_libraries ( lama_test "-Wl,--whole-archive" ama "-Wl,--no-whole-archive" ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} ${LAMA_THREAD_LIBRARY} )
    else ()
        target_link_libraries ( lama_test ama ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} ${LAMA_THREAD_LIBRARY} )
    endif ( LAMA_LIBRARY_TYPE STREQUAL "STATIC" )
endif ( WIN32 )


### Copy the README and xmltest.sh file in build directory ### 
file ( COPY README DESTINATION ${CMAKE_BINARY_DIR}/test )
file ( COPY xmltest.sh DESTINATION ${CMAKE_BINARY_DIR}/test ) 

if ( LAMA_USE_CODE_COVERAGE )
	### Copy the cca.sh and ccajen.sh file in build directory ###
	configure_file( "${CMAKE_SOURCE_DIR}/test/cca.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/cca.sh" @ONLY)
	configure_file( "${CMAKE_SOURCE_DIR}/test/ccajen.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/ccajen.sh" @ONLY)
endif ( LAMA_USE_CODE_COVERAGE )
