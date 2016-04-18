###
 # @file scai/solver/CMakeLists.txt
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
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 25.01.2016
 # @since 2.0.0
###

add_subdirectory( distributed )

scai_project ( 

	CLASSES

	HEADERS

	SOURCES

    	BiCGstabTest
    	BiCGTest
        CGTest
        CGNETest # here may be problems in convergence ?!
        CGNRTest
        CGSTest
        CommonLoggerTest
        CriterionTest
        FileLoggerTest
        GMRESTest
        InverseSolverTest
        IterationCountTest
        IterativeSolverTest
        JacobiTest
        MINRESTest
        QMRTest
        ResidualStagnationTest
        ResidualThresholdTest
        RichardsonTest
        SimpleAMGTest
        SolutionProxyTest
        SolverTest
        SolverTest1
        TFQMRTest
        TimerTest
        
   )
   
### add definitions ###
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../res/testfiles" )
   
add_executable ( solverTest ${CXX_SOURCES} )

target_link_libraries ( solverTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${PROJECT_NAME} )
    
file ( GLOB LOCAL_HEADERS "*.hpp" )

install ( FILES ${LOCAL_HEADERS} DESTINATION include/scai/solver/test )

### Copy the README and xmltest.sh file in build directory ### 
file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverage.sh" @ONLY)
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/../../scai_code_coverage_functions.sh" "${CMAKE_CURRENT_BINARY_DIR}/" @ONLY)
endif ( USE_CODE_COVERAGE )
