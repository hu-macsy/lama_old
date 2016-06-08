###
 # @file solver/test/distributed/CMakeLists.txt
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
 # @endlicense
 #
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 25.01.2016
###

scai_project ( 

	CLASSES

	HEADERS

	SOURCES
        SolverDistTest
        
        P_CGTest
        P_JacobiTest
   )
   
### set include directories ###

include_directories ( ${SCAI_BOOST_INCLUDE_DIR} )

### add definitions ###

add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../../res/testfiles" )

add_executable ( solverDistTest ${CXX_SOURCES} )

target_link_libraries ( solverDistTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${PROJECT_NAME} )
    
file ( GLOB LOCAL_HEADERS "*.hpp" )

install ( FILES ${LOCAL_HEADERS} DESTINATION include/scai/solver/test/distributed )    

### Copy the xmltest.sh file in build directory ### 
#file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )  
