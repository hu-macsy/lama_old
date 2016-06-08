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
 # @brief test/distributed/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
###

### set include directories ### 
#include_directories ( . .. )

### Add classes and header files ###

if ( METIS_FOUND OR PARMETIS_FOUND )
    set ( MetisDistributionTest MetisDistributionTest )
else( METIS_FOUND OR PARMETIS_FOUND )
    set ( MetisDistributionTest )
endif ( METIS_FOUND OR PARMETIS_FOUND )

scai_project (

   CLASSES

        P_SparseMatrixTest

   SOURCES

        P_COOSparseMatrixTest
        P_CSRSparseMatrixTest
        P_DenseMatrixTest
        P_DIASparseMatrixTest
        P_ELLSparseMatrixTest
        P_JDSSparseMatrixTest
        P_MatrixCreatorTest
        P_MatrixTest
        P_VectorTest
        P_NormTest

        ${MetisDistributionTest}

        ../SparseMatrixHelper
)
    
# individual preprocessor definition 

#set_property( SOURCE ${CXX_SOURCES} 
#              PROPERTY COMPILE_DEFINITIONS LAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../../res/testfiles" 
#            )
add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../../res/testfiles" )

### add definitions ###
if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### set link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif( WIN32 )

### add executable ###
add_executable ( lamaDistTest lama_dist_test.cpp ${CXX_SOURCES} )

### set target link libraries ###
if ( WIN32 )
    target_link_libraries ( lamaDistTest ${PROJECT_NAME} )
else ( WIN32 )
    target_link_libraries ( lamaDistTest ${PROJECT_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
endif ( WIN32 )
