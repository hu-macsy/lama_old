###
 # @file scai/hmemo/CMakeLists.txt
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
 # @brief scai/hmemo/CMakeLists.txt
 # @author Thomas Brandes
 # @date 30.06.2015
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project( scai_hmemo CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

## Define all search packages
include ( InternalDependencies )
include ( ExternalDependencies )

# Find required packages

foreach    ( module ${SCAI_HMEMO_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_HMEMO_EXTERNAL_DEPS} )

foreach    ( PACKAGE_TO_FIND ${SCAI_HMEMO_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

include ( Functions/checkValue )
include ( Functions/setAndCheckCache )

if    ( NOT SCAI_COMPLETE_BUILD )

	### Include Compiler Flags
	include ( SetCPPFlags )
	
else  ( NOT SCAI_COMPLETE_BUILD )

	add_definitions( ${ADDITIONAL_WARNING_FLAGS} )

	set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS}")
	set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS_RELEASE} " )
	set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
	set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )

endif ( NOT SCAI_COMPLETE_BUILD )

if ( NOT CXX_SUPPORTS_C11 )
    find_package( Boost REQUIRED )
    include_directories( ${Boost_INCLUDE_DIR} )
endif ( NOT CXX_SUPPORTS_C11 )

## Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_HMEMO_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_HMEMO_INTERNAL_DEPS} )

### compile flag for logging level mandatory

add_definitions( -D${SCAI_LOGGING_FLAG} )

if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

## add subdirectories
add_subdirectory ( examples )
add_subdirectory ( test )

## build library
file ( GLOB CPP_SOURCES "*.cpp" )
file ( GLOB HPP_HEADERS "*.hpp" )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
    include_directories ( ${CUDA_TOOLKIT_INCLUDE} )
endif ( CUDA_FOUND AND USE_CUDA )

include_directories ( . )

### add library and put path to search path ###

add_library ( scai_hmemo SHARED ${CPP_SOURCES} ${CUDA_SOURCES} )
target_link_libraries ( scai_hmemo ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY}
                               ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} )

foreach    ( module ${SCAI_HMEMO_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( scai_hmemo ${${upper_module}_LIBRARY} )
endforeach    ( module ${SCAI_HMEMO_INTERNAL_DEPS} )

### install ###

install ( TARGETS scai_hmemo DESTINATION lib )
install ( FILES ${HPP_HEADERS} DESTINATION include/scai/hmemo )

## install library.hpp
install ( FILES ../hmemo.hpp DESTINATION include/scai )

## add custum target distClean, but only if this project is not included

if ( CMAKE_SOURCE_DIR EQUAL CMAKE_CURRENT_SOURCE_DIR )
    add_custom_target ( distclean )
endif ()
