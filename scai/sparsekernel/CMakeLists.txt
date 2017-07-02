###
 # @file scai/sparsekernel/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief scai/sparsekernel/CMakeLists.txt
 # @author Eric Schricker
 # @date 18.02.2016
###

## Define name of the SCAI module project

set ( MODULE_NAME sparsekernel )

set ( INTERNAL_DEPS common logging tracing tasking hmemo kregistry blaskernel utilskernel ) 
set ( EXTERNAL_DEPS OpenMP CUDA MIC SCAI_BLAS ) 

## Define all library names with the (global) prefix SCAI_LIBRARY_PREFIX

set ( PROJECT_NAME "${SCAI_LIBRARY_PREFIX}${MODULE_NAME}" )

### includes ###

## important Variables for project
include ( SetProjectVariables)
include ( Dependencies/minimalSupportedVersions )

### find all external packages via the provided wrappers of CMake modules

foreach ( module ${EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ()

### finalize build flags
include ( ConcludeFlags )

### add definitions
### add definitions

include ( Settings/logging )
include ( Settings/tracing )

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_LOGGING_LEVEL_${SCAI_LOGGING_LEVEL} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )
add_definitions ( -DSCAI_TRACING_${SCAI_TRACING} )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

if    ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

## Project System-Include
include_directories ( ../.. ${CMAKE_BINARY_DIR}/include )

### add includes for external packages

foreach ( module ${EXTERNAL_DEPS} )
    string ( TOUPPER ${module} upper_module )
    include_directories( ${SCAI_${upper_module}_INCLUDE_DIR} )
endforeach ()


### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...
scai_project (

    CLASSES                  # .cpp, .hpp


    HEADERS                  # .hpp only

        COOKernelTrait
        CSRKernelTrait
        DenseKernelTrait
        DIAKernelTrait
        ELLKernelTrait
        JDSKernelTrait
        StencilKernelTrait
        
)

## add subdirectories

if    ( CUDA_FOUND AND USE_CUDA )
    add_subdirectory ( cuda )
    cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
    set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
endif ( CUDA_FOUND AND USE_CUDA )

if    ( BUILD_EXAMPLES )
add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

add_subdirectory ( external )

if ( USE_MIC )
    add_subdirectory ( mic )
endif ( USE_MIC )

add_subdirectory ( openmp )

if    ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### set link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

### add library ###

scai_library()

if ( USE_MIC )
    set ( MKLROOT $ENV{MKLROOT} )
    set ( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} "-offload-option,mic,compiler,\"-L${MKLROOT}/lib/mic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group\"" )
    #ToDo: set_target_properties not working, any better way?
endif ( USE_MIC )

### install ###

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION include/scai/${MODULE_NAME} )

## install library.hpp
install ( FILES ../${PROJECT_SURNAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

include ( CustomCommands/SphinxDoc )
