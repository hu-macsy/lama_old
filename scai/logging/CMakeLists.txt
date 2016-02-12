###
 # @file scai/logging/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief CMake configuration file for Logging library
 # @author Thomas Brandes
 # @date 09.06.2015
 # @since 2.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( scai_logging CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### includes ###

## important Variables for project
include ( VersionDefinition )
include ( SetProjectVariables)
include ( InternalDependencies )
include ( ExternalDependencies )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( CompilerVersion )
    include ( SCAIAssert )
    include ( SetCPPFlags )
else  ( NOT SCAI_COMPLETE_BUILD )
    include ( Functions/CheckValueAtCompleteBuild )
    CheckValueAtCompleteBuild ( ${PROJECT_SURNAME} )
endif ( NOT SCAI_COMPLETE_BUILD )

## C++11 or Boost
include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

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

### Project System-Include
include_directories ( ../.. )
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
add_subdirectory ( examples )
## begin subdir test
# prepare test
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/test/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/test/make.inc" )
file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/test DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
# add test
enable_testing()
add_test ( logging_test test/test.sh )
## end subdir test

### add library ###
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${${UPPER_PROJECT_NAME}_VERSION} )

## add link libraries
addInternalAndExternalLinkLibraries()

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

## install library.hpp
install ( FILES ../${PROJECT_SURNAME}.hpp DESTINATION include/scai )

## add custum target distClean
include ( DistClean )

### summary ###
if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/${PROJECT_SURNAME} )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
