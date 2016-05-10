###
 # @file CMakeLists.txt
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
 # @brief CMake configuration file for scai_blaskernel tests
 # @author Eric Schricker
 # @date 18.11.2015
###

set ( CXX_SOURCES
      BLASKernelTest
      BLAS1Test
      BLAS2Test
      BLAS3Test
      LAPACKTest
) 

include_directories ( ${SCAI_BOOST_INCLUDE_DIR} )

add_executable ( blaskernelTest ${CXX_SOURCES} )

target_link_libraries ( blaskernelTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${PROJECT_NAME} )

### Copy the xmltest.sh file in build directory ### 

file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverageFirst.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverageFirst.sh" @ONLY)
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverageSecond.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverageSecond.sh" @ONLY)
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverageThird.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverageThird.sh" @ONLY)
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/../../scai_code_coverage_functions.sh" "${CMAKE_CURRENT_BINARY_DIR}/" @ONLY)
endif ( USE_CODE_COVERAGE )
