###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
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
      BinaryOpTest
      ConstantsTest
      ContextTypeTest
      ExceptionTest  
      FactoryTest  
      Factory1Test  
      FactoryTest
      GridTest
      LibModuleTest
      MathTest
      ptrTest
      PrintableTest
      ThreadTest  
      TypeListTest
      TypeLoopTest
      TypeTraitsTest  
      UnaryOpTest
      ConstantsTest
      ScalarTypeTest  
      SettingsTest  
      WalltimeTest  
)

if    ( ${USE_COMPLEX} )
    set ( CXX_SOURCES ${CXX_SOURCES} ComplexTest TypeConversionTest )
endif ( ${USE_COMPLEX} )

include_directories ( ${BOOST_INCLUDE_DIR} )

add_executable ( commonTest ${CXX_SOURCES} )

target_link_libraries ( commonTest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${PROJECT_NAME} )

if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA )

file ( GLOB LOCAL_HEADERS "*.hpp" )

install ( FILES ${LOCAL_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR}/test )

### Copy the xmltest.sh file in build directory ### 
file ( COPY xmltest.sh DESTINATION ${CMAKE_CURRENT_BINARY_DIR} FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} ) 

if ( USE_CODE_COVERAGE )
    configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/code_coverage.sh" @ONLY)
    file ( COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/code_coverage.sh 
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
    file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/../../scai_code_coverage_functions.sh 
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
endif ( USE_CODE_COVERAGE )
