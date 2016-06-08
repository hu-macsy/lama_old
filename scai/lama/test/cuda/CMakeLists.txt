###
 # @file CMakeLists.txt
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief test/cuda/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
###

scai_project (

  CLASSES

        CUDATestContext

  SOURCES

#       CUDA_ContextTest
#       CUDA_HostContextTest
        CUDA_MatrixStorageTest
#       CUDA_VectorTest
        CUDA_KernelRegistryTest

)

# Due to additional include files an own executable is defined for CUDA tests

if ( ${SCAI_ASSERT_LEVEL} STREQUAL "OFF" )
   message ( STATUS "WARNING: tests will not check for exceptions (assertions are switched off)" )
endif ( ${SCAI_ASSERT_LEVEL} STREQUAL "OFF" )

### set include directories ###
include_directories ( ${SCAI_CUDA_INCLUDE_DIR} )

### set link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif( WIN32 )

### add executable ###
add_executable ( lamaCUDATest lama_cuda_test.cpp ${CXX_SOURCES} )

### set target link libraries ###
if ( WIN32 )
    target_link_libraries ( lamaCUDATest ${PROJECT_NAME} )
else ( WIN32 )
    target_link_libraries ( lamaCUDATest ${PROJECT_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ${Boost_REGEX_LIBRARY} )
endif ( WIN32 )
