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
 # @brief CMakeLists of common/test
 # @author Lauretta Schubert
 # @date 09.09.2015
###

set ( CXX_SOURCES
	    atomicOpenMPTest
	    bindFunctionTest
      commonTest
      ConstantsTest
      ContextTypeTest
      FactoryTest  
      Factory1Test  
      FactoryTest
      LibModuleTest
      MathTest
      ptrTest
      PrintableTest
      ThreadTest  
      TypeListTest
      TypeLoopTest
      TypeTraitsTest  
      ScalarTypeTest  
      SettingsTest  
      WalltimeTest  
)

if    ( ${USE_COMPLEX} )
    set ( CXX_SOURCES ${CXX_SOURCES} ComplexTest TypeConversionTest )
endif ( ${USE_COMPLEX} )

include_directories ( ${SCAI_BOOST_INCLUDE_DIR} )

add_executable ( commonTest ${CXX_SOURCES} )

target_link_libraries ( commonTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${PROJECT_NAME} )

if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA )

file ( GLOB LOCAL_HEADERS "*.hpp" )

install ( FILES ${LOCAL_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR}/test )

### Copy the xmltest.sh file in build directory ### 
file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}/code_coverage.sh" @ONLY)
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/../../scai_code_coverage_functions.sh" "${CMAKE_CURRENT_BINARY_DIR}/" @ONLY)
endif ( USE_CODE_COVERAGE )
