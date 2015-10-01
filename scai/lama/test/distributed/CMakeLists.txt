###
 # @file CMakeLists.txt
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
 # @brief test/distributed/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
 # @since 1.0.0
###

### set this dir as lama source dir ###
scai_project_start ()

### set include directories ### 
#include_directories ( . .. )

### Add classes and header files ###

if ( METIS_FOUND OR PARMETIS_FOUND )
    set ( MetisDistributionTest MetisDistributionTest )
else( METIS_FOUND OR PARMETIS_FOUND )
    set ( MetisDistributionTest )
endif ( METIS_FOUND OR PARMETIS_FOUND )

scai_project_add(

   CLASSES

        DistributionTest
        P_SparseMatrixTest

   SOURCES

        BlockDistributionTest
        CyclicDistributionTest
        GenBlockDistributionTest
        GeneralDistributionTest
        NoDistributionTest
        P_CGTest
        P_COOSparseMatrixTest
        P_CSRSparseMatrixTest
        P_DenseMatrixTest
        P_DIASparseMatrixTest
        P_ELLSparseMatrixTest
        P_InverseSolverTest
        P_JacobiTest  
        P_JDSSparseMatrixTest
        P_MatrixCreatorTest
        P_MatrixStorageTest
        P_MatrixTest
        P_MPICommunicatorTest
        P_RedistributorTest
        P_SORTest
        P_SpecializedJacobiTest
        P_VectorTest
        P_NormTest

        ${MetisDistributionTest}

        #TODO: REMOVE!
        ../Configuration
        ../SameMatrixHelper
        ../SparseMatrixHelper
        ../CommunicatorTest
)
    
# individual preprocessor definition 

set_property( SOURCE ../Configuration.cpp 
              PROPERTY COMPILE_DEFINITIONS LAMA_TESTFILE_PATH="${CMAKE_SOURCE_DIR}/../res/testfiles" 
            )

### add definitions ###
if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### set link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif( WIN32 )

### add executable ###
add_executable ( lama_dist_test lama_dist_test.cpp ${CXX_SOURCES} )

### set target link libraries ###
if ( WIN32 )
    target_link_libraries ( lama_dist_test ama )
else ( WIN32 )
    target_link_libraries ( lama_dist_test ama ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} )
endif ( WIN32 )
