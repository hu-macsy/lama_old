###
 # @file solver/test/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @endlicense
 #
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 25.01.2016
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
