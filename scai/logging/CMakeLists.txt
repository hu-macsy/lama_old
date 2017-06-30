###
 # @file scai/logging/CMakeLists.txt
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
 # @brief CMake configuration file for Logging library
 # @author Thomas Brandes
 # @date 09.06.2015
###

### includes ###

project( scai_logging )

## important Variables for project
set ( SCAI_LOGGING_VERSION    1.0.0 )
include ( SetProjectVariables)
set ( SCAI_LOGGING_EXTERNAL_DEPS                                                          ) # 2
include ( Dependencies/minimalSupportedVersions )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
include ( SetBuildFlags )
include ( VersionCheck/Compiler )
include ( Compiler/SetCPPFlags )

## Boost for Test
if ( BUILD_TEST )
    include( Package/BoostTest )
endif ( BUILD_TEST )

### include all search packages
addInternalAndExternalPackages()

## finalize build flags
include ( ConcludeFlags )

### add definitions
add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )
# logging level for logging library
# GenLogger uses also logging during initialization, here we make a fixed setting
add_definitions ( -DSCAI_LOG_LEVEL_DEBUG )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include ( also configured files )
set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )
addInternalAndExternalIncludes()

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        GenLoggerCreator
        Level
        Logger
        LoggerProvider
        GenLogger
        SourceLocation
 
    HEADERS                  # .hpp only

        AbstractLoggerCreator
)

## add subdirectories
add_subdirectory ( doc )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

add_subdirectory ( macros )
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

set ( TEST_FILES Makefile README complexLogging.cpp simpleLogging.cpp test.sh )
foreach    ( TEST_FILE ${TEST_FILES} )
    file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/test/${TEST_FILE}
           DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/test
           FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
endforeach ( TEST_FILE ${TEST_FILES} )

# add test

enable_testing()
add_test ( logging_test test/test.sh )
## end subdir test

### set link directories ###
if ( WIN32 )
  link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add library ###
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${${UPPER_PROJECT_NAME}_VERSION} )

## ToDo: add link libraries
# addInternalAndExternalLinkLibraries()
target_link_libraries( ${PROJECT_NAME} scai_common )

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

## install library.hpp
install ( FILES ../${PROJECT_SURNAME}.hpp DESTINATION include/scai )

## add custum target distClean
include ( CustomCommands/DistClean )
