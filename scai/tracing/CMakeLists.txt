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
 # @brief tracing/CMakeLists.txt
 # @author Thomas Brandes
 # @date 10.06.2015
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project( scai_tracing CXX C )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

## Define all search packages
include ( InternalDependencies )
include ( SearchPackages/Tracing ) 

# Find required packages
foreach    ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} QUIET REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
foreach    ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )

include ( Functions/checkValue )
include ( Functions/setAndCheckCache )

include ( Package/CUDA )

### Include Compiler Flags
include ( SetBuildFlags )
include ( SetCPPFlags )

## Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_TRACING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_TRACING_INTERNAL_DEPS} )

### Mandatory: compile flag for logging level

add_definitions ( -D${SCAI_LOGGING_FLAG} )

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

## add subdirectory
add_subdirectory ( examples )

## build library
file ( GLOB CPP_SOURCES "*.cpp" )
file ( GLOB HPP_HEADERS "*.hpp" )

### add library and put path to search path ###

add_library ( scai_tracing SHARED ${CPP_SOURCES} )

foreach    ( module ${SCAI_TRACING_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( scai_tracing ${${upper_module}_LIBRARY} )
endforeach    ( module ${SCAI_TRACING_INTERNAL_DEPS} )

### install ###

install ( TARGETS scai_tracing DESTINATION lib/ )
install ( FILES ${HPP_HEADERS} DESTINATION include/scai/tracing )

if ( CMAKE_VERSION VERSION_GREATER 2.8.7 ) 
	add_subdirectory( gui )
else ()
	message( STATUS "GUI will not be build" )
endif( )

## install library.hpp
install ( FILES ../tracing.hpp DESTINATION include/scai )

## add custum target distClean, but only if this project is not included

if ( CMAKE_SOURCE_DIR EQUAL CMAKE_CURRENT_SOURCE_DIR )
    add_custom_target ( distclean )
endif ()
