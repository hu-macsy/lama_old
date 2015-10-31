###
 # @file scai/common/CMakeLists.txt
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
 # @brief CMake configuration file for common library
 # @author Thomas Brandes
 # @date 09.06.2015
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( scai_common CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### Include Compiler Flags
if    ( NOT SCAI_COMPLETE_BUILD )
	include ( SetBuildFlags )
	include ( SCAIAssert )
	include ( SetCPPFlags )
else  ( NOT SCAI_COMPLETE_BUILD )
	include ( Functions/CheckValueAtCompleteBuild )
	CheckValueAtCompleteBuild ( common )
endif ( NOT SCAI_COMPLETE_BUILD )

include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 )

### include all search packages
# no internal dependencies
include ( ExternalDependencies )

foreach    ( module ${SCAI_COMMON_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_COMMON_EXTERNAL_DEPS} )

### finalize build flags
include ( ConcludeFlags )

### add definitions

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include
include_directories ( ../.. )

### build library

## add subdirectories
add_subdirectory ( examples )

## add source files

include ( Functions/scaiProject )

scai_project (

    CLASSES

        Printable
        ScalarType
        SCAITypes
        ContextType
        Settings
        Thread
        Walltime

    HEADERS

        Assert
        bind     
        Complex
        config
        Constants
        Factory
        Factory1
        function
        OpenMP
        NonCopyable
        shared_ptr
        TypeTraits
        unique_ptr
        weak_ptr
)

## add CUDA if supported
if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
    include_directories ( ${CUDA_TOOLKIT_INCLUDE} )
endif ( CUDA_FOUND AND USE_CUDA )

# subdirectory macros
add_subdirectory ( macros )

## add exception dir
add_subdirectory ( exception )

add_subdirectory ( test )

## add library
add_library ( scai_common ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

## add link libraries
foreach    ( module ${SCAI_COMMON_EXTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( scai_common ${SCAI_${upper_module}_LIBRARIES} )
endforeach ( module ${SCAI_COMMON_EXTERNAL_DEPS} )

### install ###
install ( TARGETS scai_common DESTINATION lib )

include ( Functions/relativeInstall )
relative_install ( FILES ${CXX_HEADERS} DESTINATION include/scai/common )

## install library.hpp
install ( FILES ../common.hpp DESTINATION include/scai )

## add custum target distClean
include ( DistClean )

if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
	include ( Summaries/common )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
