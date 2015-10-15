###
 # @file scai/kernel/CMakeLists.txt
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
 # @brief scai/kernel/CMakeLists.txt
 # @author Thomas Brandes
 # @date 30.06.2015
 # @since 2.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project( scai_kernel CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### Include Compiler Flags
if    ( NOT SCAI_COMPLETE_BUILD )
	include ( SetBuildFlags )
	include ( SCAIAssert )
	include ( SetCPPFlags )
else  ( NOT SCAI_COMPLETE_BUILD )	
	include ( Functions/CheckValueAtCompleteBuild )
	CheckValueAtCompleteBuild ( kernel )
endif ( NOT SCAI_COMPLETE_BUILD )

### include all search packages
include ( InternalDependencies )
include ( ExternalDependencies )

foreach    ( module ${SCAI_KERNEL_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_KERNEL_EXTERNAL_DEPS} )

foreach    ( PACKAGE_TO_FIND ${SCAI_KERNEL_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

## finalize build flags
include ( ConcludeFlags )

### add definitions

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )

add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
add_definitions ( -D${SCAI_TRACING_FLAG} )

add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

## Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_KERNEL_INTERNAL_DEPS} ${SCAI_KERNEL_EXTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_KERNEL_INTERNAL_DEPS} ${SCAI_KERNEL_EXTERNAL_DEPS} )

### build library

include ( Functions/scaiProject )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        ContextFunction
        KernelRegistry
)

## add CUDA if supported
if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
    include_directories ( ${CUDA_TOOLKIT_INCLUDE} )
endif ( CUDA_FOUND AND USE_CUDA )

## add subdirectories
# add_subdirectory ( examples )
# add_subdirectory ( test )

## add library
add_library ( scai_kernel ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

## add link libraries
list ( REVERSE SCAI_KERNEL_INTERNAL_DEPS )
foreach    ( module ${SCAI_KERNEL_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
	set ( KERNEL_LINK_LIBRARIES ${KERNEL_LINK_LIBRARIES} ${${upper_module}_LIBRARY} )
endforeach ( module ${SCAI_KERNEL_INTERNAL_DEPS} )
target_link_libraries ( scai_kernel ${SCAI_START_LINK_LIBRARIES} ${KERNEL_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES} )

foreach    ( module ${SCAI_KERNEL_EXTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( scai_kernel ${SCAI_${upper_module}_LIBRARIES} )
endforeach ( module ${SCAI_KERNEL_EXTERNAL_DEPS} )

### install

install ( TARGETS scai_kernel DESTINATION lib )

include( Functions/relativeInstall )
relative_install( FILES ${CXX_HEADERS} DESTINATION include/scai/kernel )

## install library.hpp
install ( FILES ../kernel.hpp DESTINATION include/scai )

## add custum target distClean
include ( DistClean )

