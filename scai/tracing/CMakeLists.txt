###
 # @file scai/tracing/CMakeLists.txt
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
 # @brief scai/tracing/CMakeLists.txt
 # @author Thomas Brandes
 # @date 10.06.2015
###

## Define SCAI module project by name and internal/external dependencies

scai_module(

    MODULE_NAME   tracing
    INTERNAL_DEPS common logging
    EXTERNAL_DEPS 
)

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        CallTreeTable
        FileTable
        RegionTable
        TraceConfig
        TraceData
        TraceRegionRecord
        VTInterface

    HEADERS                  # .hpp only

        CallStack
        Counters
        RegionEntry
)

## add subdirectories

scai_subdirectories ( EXAMPLES )

## begin subdir test
# prepare test
if ( USE_CODE_COVERAGE )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/test/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/test/code_coverage.sh" @ONLY)
    file ( COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/test/code_coverage.sh 
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
    file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/../scai_code_coverage_functions.sh 
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
endif ( USE_CODE_COVERAGE )

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/test/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/test/make.inc" )

# The shell script test.sh uses NUM_THREADS_TEST to verify tracing for threads

if ( USE_OPENMP )
   set ( NUM_THREADS_TEST 4 )
else ( USE_OPENMP )
   set ( NUM_THREADS_TEST 1 )
endif ( USE_OPENMP )

# Now configure test.sh and set correct permission

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/test/test.sh.in"
                 "${CMAKE_CURRENT_BINARY_DIR}/test/conf/test.sh" @ONLY )

file ( COPY "${CMAKE_CURRENT_BINARY_DIR}/test/conf/test.sh" 
       DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/test"
       FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )

set ( TEST_FILES Makefile README simpleTracing.cpp )

foreach    ( TEST_FILE ${TEST_FILES} )
    file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_FILE}
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
endforeach ( TEST_FILE ${TEST_FILES} )

# add test
enable_testing()
add_test ( tracing_test test/test.sh )
## end subdir test

add_subdirectory( gui )
add_subdirectory ( macros )

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ()

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install library.hpp
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )
