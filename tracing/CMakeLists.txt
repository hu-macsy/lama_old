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
project( tracing CXX C )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

## Define all search packages
include ( SearchPackages/Tracing ) 

# Find required packages
foreach    ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} REQUIRED )
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
include_directories ( . )

set ( DEPS COMMON LOGGING )
foreach    ( module ${DEPS} )
    include_directories( ${${module}_INCLUDE_DIR} )
endforeach ( module ${DEPS} )

### add definitions ###

if    ( LAMA_LOG_LEVEL )
    add_definitions ( -DLAMA_LOG_LEVEL_${LAMA_LOG_LEVEL} )
else  ( LAMA_LOG_LEVEL )
    add_definitions( -DLAMA_LOG_LEVEL_DEBUG )
endif ( LAMA_LOG_LEVEL ) 

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

if ( FALSE )

### VampirTrace is not easy to add for tracing, is in development

    find_package ( VampirTrace )

    message( STATUS "VAMPIRTRACE found = ${VAMPIRTRACE_FOUND}" )
    message( STATUS "MPI found = ${MPI_FOUND}" )

    if ( ${VAMPIRTRACE_FOUND} )
        message( STATUS "VAMPIRTRACE found, lib = ${VAMPIRTRACE_LIBRARIES}" )
        add_definitions( -DUSE_VAMPIRTRACE )
        include_directories( ${VAMPIRTRACE_INCLUDE_DIR} )
        target_link_libraries( tracing ${VAMPIRTRACE_LIBRARIES} ${MPI_LIBRARIES} )
    endif ()

endif ()

add_subdirectory ( tracing )
add_subdirectory ( examples )

if ( CMAKE_VERSION VERSION_GREATER 2.8.7 ) 
	add_subdirectory( gui )
else ()
	message( STATUS "GUI will not be build" )
endif( )

## install library.hpp
install ( FILES tracing.hpp DESTINATION include )

## add custum target distClean
include ( DistClean )