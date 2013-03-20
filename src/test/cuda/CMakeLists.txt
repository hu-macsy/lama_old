###################################
#
#  Executable: lama_cuda_btests
#
###################################

# Due to additional include files an own executable is defined for CUDA tests

if ( ${LAMA_ASSERT_LEVEL} STREQUAL "OFF" )
   message ( STATUS "WARNING: tests will not check for exceptions (assertions are switched off)" )
endif ( ${LAMA_ASSERT_LEVEL} STREQUAL "OFF" )

### set include directories ### 
    
include_directories ( ${CUDA_INCLUDE_DIRS} )


### set sources ### 

set ( CUDA_TESTSOURCES
           CUDA_ContextTest.cpp
           CUDA_HostContextTest.cpp
           CUDA_LAMAArrayTest.cpp
           CUDA_MatrixStorageTest.cpp
           CUDA_VectorTest.cpp #move distributed functionality to distribute
           CUDAContext.cpp
    )
    
set ( CUDA_TESTSOURCES ${CUDA_TESTSOURCES}
      distributed/CUDA_P_SparseMatrixTest.cpp
    )  

### set link directories ###

if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif( WIN32 )

### add executable ###

add_executable ( lama_cuda_test lama_cuda_test.cpp ${CUDA_TESTSOURCES} )


### set target link libraries ###

if ( WIN32 )
    target_link_libraries ( lama_cuda_test lama )
else ( WIN32 )
    target_link_libraries ( lama_cuda_test lama ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
endif ( WIN32 )
