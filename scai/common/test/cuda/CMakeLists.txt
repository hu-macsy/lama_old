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
 # @endlicense
 #
 # @brief CMake configuration file for cuda memory tests
 # @author Thomas Brandes
 # @date 16.06.2015
###

scai_project (

    SOURCES

        commonCUDATest
        CUDAAccessTest
        CUDADeviceTest
        CUDASettingsTest
 
    CUDA_SOURCES

        CUDAKernel
)

include_directories ( ${SCAI_BOOST_INCLUDE_DIR} )

cuda_add_executable ( commonCUDATest ${CUDA_SOURCES} ${CXX_SOURCES} )

target_link_libraries ( commonCUDATest ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} scai_common )
