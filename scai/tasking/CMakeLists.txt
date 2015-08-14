###
 # @file CMakeLists.txt
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
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project( tasking CXX C )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../CMake/Modules" )

## Define all search packages
include ( SearchPackages/Tasking ) 

# Find required packages
foreach    ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
foreach    ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )

### Include Compiler Flags
include ( SetBuildFlags )
include ( SetCPPFlags )

## Project System-Include
include_directories ( ../.. )

set ( DEPS SCAI_COMMON SCAI_LOGGING SCAI_TRACING )
foreach    ( module ${DEPS} )
    include_directories( ${${module}_INCLUDE_DIR} )
endforeach ( module ${DEPS} )

### add definitions ###

## log level
if    ( SCAI_LOG_LEVEL )
    add_definitions ( -DSCAI_LOG_LEVEL_${SCAI_LOG_LEVEL} )
else  ( SCAI_LOG_LEVEL )
    add_definitions( -DSCAI_LOG_LEVEL_DEBUG )
endif ( SCAI_LOG_LEVEL )

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

add_definitions( -DSCAI_TRACE_ON )

## add subdirectories
add_subdirectory ( examples )

## build library
file ( GLOB CPP_SOURCES "*.cpp" )
file ( GLOB HPP_HEADERS "*.hpp" )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

include_directories ( . )

### add library and put path to search path ###

add_library ( scai_tasking SHARED ${CPP_SOURCES} )

foreach    ( module ${DEPS} )
    target_link_libraries ( scai_tasking ${${module}_LIBRARY} )
endforeach    ( module ${DEPS} )

### install ###

install ( TARGETS scai_tasking DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/scai/tasking )

## install library.hpp
install ( FILES ../tasking.hpp DESTINATION include/scai )

## add custum target distClean
include ( DistClean )