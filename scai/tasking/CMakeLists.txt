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
project( scai_tasking CXX )

### Set version of tasking
# TOOD: where is the best place to define this?
set ( SCAI_TASKING_VERSION_NAME "Alpamayo" )
set ( SCAI_TASKING_COPYRIGHT_YEAR 2015 )
set ( SCAI_TASKING_VERSION_MAJOR 1 )
set ( SCAI_TASKING_VERSION_MINOR 0 )
set ( SCAI_TASKING_VERSION_PATCH 1 )
set ( SCAI_TASKING_VERSION_COUNT 3 )
set ( SCAI_TASKING_VERSION ${SCAI_TASKING_VERSION_MAJOR}.${SCAI_TASKING_VERSION_MINOR}.${SCAI_TASKING_VERSION_PATCH} )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### Include Compiler Flags
if    ( NOT SCAI_COMPLETE_BUILD )
	include ( SetBuildFlags )
	include ( SCAIAssert )
	include ( SetCPPFlags )
	include ( Package/Sphinx )
else  ( NOT SCAI_COMPLETE_BUILD )
	include ( Functions/CheckValueAtCompleteBuild )
	CheckValueAtCompleteBuild ( tasking )
endif ( NOT SCAI_COMPLETE_BUILD )

include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

if    ( NOT CXX_SUPPORTS_C11 )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../boost_check.hpp.in"
                     "${CMAKE_CURRENT_BINARY_DIR}/boost_check.hpp" )
else  ( NOT CXX_SUPPORTS_C11 )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../boost_check.hpp.in_empty"
                     "${CMAKE_CURRENT_BINARY_DIR}/boost_check.hpp" )
endif ( NOT CXX_SUPPORTS_C11 )
install ( FILES ${CMAKE_CURRENT_BINARY_DIR}/boost_check.hpp DESTINATION include/scai )

### include all search packages
include ( InternalDependencies )
#no external dependencies

foreach    ( PACKAGE_TO_FIND ${SCAI_TASKING_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

## finalize build flags
include ( ConcludeFlags )

### add definitions

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )

add_definitions ( -D${SCAI_LOGGING_FLAG} )

add_definitions( -DSCAI_TRACE_ON )

add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_TASKING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_TASKING_INTERNAL_DEPS} )

### build library

## add subdirectories
add_subdirectory ( examples )

include( Functions/scaiProject )

scai_project (

    HEADERS

        NoSyncToken

    CLASSES

        SyncToken
        TaskSyncToken
        Task
        ThreadPool
)

### add library
add_library ( scai_tasking ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

## add link libraries
list ( REVERSE SCAI_TASKING_INTERNAL_DEPS )
foreach    ( module ${SCAI_TASKING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
	set ( TASKING_LINK_LIBRARIES ${TASKING_LINK_LIBRARIES} ${${upper_module}_LIBRARY} )
endforeach    ( module ${SCAI_TASKING_INTERNAL_DEPS} )
target_link_libraries ( scai_tasking ${SCAI_START_LINK_LIBRARIES} ${TASKING_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES} )

### install ###
install ( TARGETS scai_tasking DESTINATION lib )

include( Functions/relativeInstall )

relative_install ( FILES ${CXX_HEADERS} DESTINATION include/scai/tasking )

## install library.hpp
install ( FILES ../tasking.hpp DESTINATION include/scai )

### add documentation
add_subdirectory( doc )

## add custum target distClean
include ( DistClean )

if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/tasking )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
