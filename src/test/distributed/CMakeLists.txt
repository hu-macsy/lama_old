### set include directories ### 

include_directories ( . .. )

### set source files ###

set ( DIST_TESTSOURCES

        ../Configuration.cpp
        ../SameMatrixHelper.cpp
        ../SparseMatrixHelper.cpp
        ../CommunicatorTest.cpp
                        
        BlockDistributionTest.cpp
        CyclicDistributionTest.cpp
        DistributionTest.cpp
        GenBlockDistributionTest.cpp
        GeneralDistributionTest.cpp
        NoDistributionTest.cpp

        P_CGTest.cpp
        P_COOSparseMatrixTest.cpp
        P_CSRSparseMatrixTest.cpp
        P_DenseMatrixTest.cpp
        P_DIASparseMatrixTest.cpp
        P_ELLSparseMatrixTest.cpp
        P_JDSSparseMatrixTest.cpp
        P_InverseSolverTest.cpp
        P_JacobiTest.cpp
        P_JDSSparseMatrixTest.cpp
        P_MatrixCreatorTest.cpp
        P_MatrixStorageTest.cpp
        P_MatrixTest.cpp
        P_MPICommunicatorTest.cpp
        P_RedistributorTest.cpp
        #P_ShiftTest.cpp
        P_SORTest.cpp
        P_SparseMatrixTest.cpp
        P_SpecializedJacobiTest.cpp
        P_VectorTest.cpp
    )


### add definitions ###
       
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_SOURCE_DIR}/../res/testfiles" )


### set link directories ###

if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif( WIN32 )

### add executable ###

add_executable ( lama_dist_test lama_dist_test.cpp ${DIST_TESTSOURCES} )

### set target link libraries ###

if ( WIN32 )
    target_link_libraries ( lama_dist_test lama )
else ( WIN32 )
    target_link_libraries ( lama_dist_test lama ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} )
endif ( WIN32 )

### add definitions ###

if ( WIN32 )
    add_definitions ( -DLAMA_COMPILING_DLL )
endif ( WIN32 )
