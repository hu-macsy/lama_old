########################################################################
#
#  CMakeLists.txt : CMake configuration for unit tests
#
########################################################################

if ( ${LAMA_ASSERT_LEVEL} STREQUAL "OFF" )
   message ( STATUS "WARNING: tests will not check for exceptions (assertions are switched off)" )
endif ( ${LAMA_ASSERT_LEVEL} STREQUAL "OFF" )

### set include directories ### 

if ( NOT ( Boost_UNIT_TEST_FRAMEWORK_LIBRARY OR Boost_REGEX ) )
    message ( FATAL_ERROR "No Boost unit test available." )
endif ( NOT ( Boost_UNIT_TEST_FRAMEWORK_LIBRARY OR Boost_REGEX ) )

include_directories ( .. ../lama ${Boost_INCLUDE_DIRS} )


### set source files ###

set ( TESTSOURCES 
        CriterionTest.cpp
        CGTest.cpp
        CommonLoggerTest.cpp
        CommunicatorTest.cpp
        Configuration.cpp
        Context.cpp
        COOSparseMatrixTest.cpp
        COOStorageTest.cpp
        CSRSparseMatrixTest.cpp
        CSRStorageTest.cpp
        DenseMatrixTest.cpp
        DenseMatrixTest1.cpp
        DenseStorageTest.cpp
        DIASparseMatrixTest.cpp
        DerivedSparseMatrixTest.cpp
        DIAStorageTest.cpp
        ELLSparseMatrixTest.cpp
        ELLStorageTest.cpp
        ELLUtilsTest.cpp
        ExceptionTest.cpp
        #ExpressionTest.cpp
        FileLoggerTest.cpp
        GMRESTest.cpp
        HostReadAccessTest.cpp
        HostWriteAccessTest.cpp
        InverseSolverTest.cpp
        IterationCountTest.cpp
        JacobiTest.cpp
        JDSSparseMatrixTest.cpp
        JDSStorageTest.cpp
        JDSUtilsTest.cpp
        L1NormTest.cpp
        L2NormTest.cpp
        LAMAArrayTest.cpp
        #LUSolverTest.cpp
        MatrixStorageConversionTest.cpp
        MatrixStorageTest.cpp
        MatrixTest.cpp
        MaxNormTest.cpp
        MetaSolverTest.cpp
        NoCommunicatorTest.cpp
        NormTest.cpp
        OpenMPTimerTest.cpp
        ResidualStagnationTest.cpp
        ResidualThresholdTest.cpp
        SameMatrixHelper.cpp
        ScalarTest.cpp
        SimpleAMGTest.cpp
        SolutionProxyTest.cpp
        SORTest.cpp
        SparseAssemblyStorageTest.cpp
        SparseMatrixHelper.cpp
        SparseMatrixTest.cpp
        SpecializedJacobiTest.cpp
        StorageIOTest.cpp
        ThreadPoolTest.cpp
        TaskSyncTokenTest.cpp
        #TimerTest.cpp
        TracingTest.cpp
        UtilsTest.cpp
        VectorTest.cpp
    )

    
### add subdirectories ###

if ( MPI_FOUND OR ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) )
    add_subdirectory ( distributed )
endif ( MPI_FOUND OR ( ( OSHMEM_FOUND AND GASNET_FOUND) OR GPI_FOUND ) )

if ( CUDA_FOUND )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND )

if ( OPENCL_FOUND )
    add_subdirectory ( OpenCL )
endif ( OPENCL_FOUND )
 

### add definitions ###
       
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_SOURCE_DIR}/../res/testfiles" )


### add link directories ###

if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )


### add executeable ###

add_executable ( lama_test lama_test.cpp ${TESTSOURCES} )

### target link libraries ###

if ( WIN32 )
    target_link_libraries ( lama_test lama )
else ( WIN32 )
    target_link_libraries ( lama_test lama ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} )
endif ( WIN32 )

### Copy the README file in build directory ### 

file ( COPY README DESTINATION ${CMAKE_BINARY_DIR}/test )

if ( CODE_COVERAGE )
	### Copy the cca.sh file in build directory ###
	file ( COPY cca.sh DESTINATION ${CMAKE_BINARY_DIR}/test ) 
endif ( CODE_COVERAGE )
