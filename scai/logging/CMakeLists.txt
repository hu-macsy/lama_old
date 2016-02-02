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
project( scai_logging CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

include ( VersionDefinition )

### Include Compiler Flags
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( CompilerVersion )
    include ( SCAIAssert )
    include ( SetCPPFlags )
    include ( Package/Sphinx )
else  ( NOT SCAI_COMPLETE_BUILD )
	include ( Functions/CheckValueAtCompleteBuild )
	CheckValueAtCompleteBuild ( logging )
endif ( NOT SCAI_COMPLETE_BUILD )

include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

### include all search packages
include ( InternalDependencies )
# no external dependencies

foreach    ( PACKAGE_TO_FIND ${SCAI_LOGGING_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${SCAI_LOGGING_INTERNAL_DEPS} )

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

foreach    ( module ${SCAI_LOGGING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_LOGGING_INTERNAL_DEPS} )

### build library

## add subdirectories
add_subdirectory ( examples )

## set files in variables CXX_SOURCES and CXX_HEADERS

include ( Functions/scaiProject )

scai_project (

    CLASSES

        GenLoggerCreator
        Level
        Logger
        LoggerProvider
        GenLogger
        SourceLocation
 
    HEADERS

        AbstractLoggerCreator
)

### add library

add_library ( scai_logging ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

## add link libraries
list ( REVERSE SCAI_LOGGING_INTERNAL_DEPS )
foreach    ( module ${SCAI_LOGGING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
	set ( LOGGING_LINK_LIBRARIES ${LOGGING_LINK_LIBRARIES} ${${upper_module}_LIBRARY} )
endforeach ( module ${SCAI_LOGGING_INTERNAL_DEPS} )
target_link_libraries ( scai_logging ${SCAI_START_LINK_LIBRARIES} ${LOGGING_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES} )

### install ###
install ( TARGETS scai_logging DESTINATION lib )

include ( Functions/relativeInstall )
relative_install ( FILES ${CXX_HEADERS} DESTINATION include/scai/logging )

## install library.hpp
install ( FILES ../logging.hpp DESTINATION include/scai )

## prepare test
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/test/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/test/make.inc" )
file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/test DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )

## add test
enable_testing()
add_test ( logging_test test/test.sh )

### add documentation
add_subdirectory( doc )

## add custum target distClean
include ( DistClean )

if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/logging )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
