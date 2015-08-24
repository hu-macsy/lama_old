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

if    ( NOT SCAI_COMPLETE_BUILD )
	
	### Include Compiler Flags
	include ( SetBuildFlags )
	include ( SetCPPFlags )
	
endif ( NOT SCAI_COMPLETE_BUILD )

if ( NOT CXX_SUPPORTS_C11 )
    find_package( Boost REQUIRED )
    include_directories( ${Boost_INCLUDE_DIR} )
endif ( NOT CXX_SUPPORTS_C11 )

## Define all search packages
include ( InternalDependencies )
#no external dependencies

# Find required packages

foreach    ( PACKAGE_TO_FIND ${SCAI_LOGGING_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${SCAI_LOGGING_INTERNAL_DEPS} )

## Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_LOGGING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_LOGGING_INTERNAL_DEPS} )

### add definitions ###

# GenLogger uses also logging during initialization, here we make a fixed setting

add_definitions ( -DSCAI_LOG_LEVEL_DEBUG )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

add_subdirectory ( examples )

## build library
file ( GLOB CPP_SOURCES "*.cpp" )
file ( GLOB HPP_HEADERS "*.hpp" )

#### concluding all defined compiler flags to CMAKE_..._FLAGS ####

## add variables to cache with new names so they can be modified by the user via CCMAKE

set ( ADDITIONAL_CXX_FLAGS "${LAMA_CXX_FLAGS}" CACHE STRING "additional flags for cxx compile and link" )
set ( ADDITIONAL_WARNING_FLAGS "${LAMA_WARNING_FLAGS}" CACHE STRING "compilation flags concerning warnings" )
set ( ADDITIONAL_CXX_FLAGS_RELEASE "${LAMA_CXX_FLAGS_RELEASE}" CACHE STRING "addtional cxx compiler flags for release optimizations" )
set ( ADDITIONAL_LINKER_FLAGS "${LAMA_LINKER_FLAGS}" CACHE STRING "additional linker flags" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS ADDITIONAL_WARNING_FLAGS ADDITIONAL_CXX_FLAGS_RELEASE ADDITIONAL_LINKER_FLAGS )

###  Code coverage with gcov/lcov
set ( USE_CODE_COVERAGE FALSE CACHE BOOL "Enable / Disable use of Code Coverage" )

if    ( USE_CODE_COVERAGE )
    set ( COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
endif ( USE_CODE_COVERAGE )

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS}")
set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS_RELEASE} " )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )

### add definitions ###

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )

### add library and put path to search path ###

add_library ( scai_logging SHARED ${CPP_SOURCES} )

foreach    ( module ${SCAI_LOGGING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( scai_logging ${${upper_module}_LIBRARY} )
endforeach    ( module ${SCAI_LOGGING_INTERNAL_DEPS} )

### install ###

install ( TARGETS scai_logging DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/scai/logging )

## install library.hpp
install ( FILES ../logging.hpp DESTINATION include/scai )

## add custum target distClean, but only if this project is not included

if ( CMAKE_SOURCE_DIR EQUAL CMAKE_CURRENT_SOURCE_DIR )
    add_custom_target ( distclean )
endif ()
