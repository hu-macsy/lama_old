###
 # @file scai/tasking/CMakeLists.txt
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
 # @brief CMake configuration file for tasking library
 # @author Thomas Brandes
 # @date 10.06.2015
 # @since 2.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( scai_tasking CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### includes ###

## important Variables for project
include ( VersionDefinition )
include ( SetProjectVariables)
include ( Dependencies/internal )
include ( Dependencies/external )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( Compiler/VersionCheck )
    include ( Compiler/SetCPPFlags )
else  ( NOT SCAI_COMPLETE_BUILD )
    include ( Functions/CheckValueAtCompleteBuild )
    CheckValueAtCompleteBuild ( ${PROJECT_SURNAME} )
endif ( NOT SCAI_COMPLETE_BUILD )

## Use Boost if BUILD_TEST is enabled
if ( BUILD_TEST )
    include( Package/Boost )
endif ( BUILD_TEST )

### include all search packages
addInternalAndExternalPackages()

## finalize build flags
include ( ConcludeFlags )
include ( SetCacheVariableTexts )

### add definitions

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions( -DSCAI_TRACE_ON )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include ( also configured files )
set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )
addInternalAndExternalIncludes()

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    HEADERS                  # .hpp only

        NoSyncToken

    CLASSES                  # .cpp, .hpp

        SyncToken
        TaskSyncToken
        Task
        ThreadPool
)

## add subdirectories

if ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA )

if    ( USE_MIC )
    add_subdirectory ( mic )
endif ( USE_MIC )

add_subdirectory ( doc )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

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
include ( CustomCommands/DistClean )

### summary ###
if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/${PROJECT_SURNAME} )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
