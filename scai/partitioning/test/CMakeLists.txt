###
 # @file scai/partitioning/test/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake configuration file for SCAI module partitioning
 # @author Thomas Brandes
 # @date 04.07.2017
###

set ( TEST_SOURCES PartitioningTest partitioningTest )

if ( METIS_FOUND AND USE_METIS ) 
    set ( TEST_SOURCES ${TEST_SOURCES} MetisPartitioningTest )
endif ()

if ( PARMETIS_FOUND AND USE_PARMETIS ) 
    # set ( TEST_SOURCES ${TEST_SOURCES} ParMetisTest )
endif ()

### set this dir as lama source dir ###

scai_project ( 

    SOURCES

        ${TEST_SOURCES}

    HEADERS

        TestDistributions
    )

scai_test ( EXECUTABLE partitioningTest 
            FILES      ${CXX_SOURCES} 
            UNIT_TEST
            RUN
)   

###  copy xmltest.sh script

scai_test_scripts( SCRIPTS       xmltest.sh )

###  configure code_coverage.sh (if USE_CODE_COVERAGE ) in build direcotry

scai_test_code_coverage( ${USE_CODE_COVERAGE} )
