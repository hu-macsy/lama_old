###
 # @file scai/lama/test/storage/CMakeLists.txt
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
 # @brief Build configuration for test of LAMA storage classes
 # @author Thomas Brandes
 # @date 03.07.2017
###

scai_project ( 

    CLASSES
    
    HEADERS

        Storages
        TestStorages
        StorageTemplateTests
        
    SOURCES

        AllStorageTest
        TypedStorageTest
        StorageConstructorTest

        COOUtilsTest

        COOStorageTest
        CSRStorageTest
        ELLStorageTest
        DIAStorageTest
        JDSStorageTest
        StencilStorageTest
        DenseStorageTest

        IOStreamTest
        FileIOTest
        MMIOTest

        lamaStorageTest
    )

### add definitions ###

add_definitions ( -DLAMA_TESTFILE_PATH="${CMAKE_CURRENT_SOURCE_DIR}/../../../res/testfiles" )

scai_test ( EXECUTABLE lamaStorageTest 
            FILES      ${CXX_SOURCES} 
            RUN 
            UNIT_TEST )
