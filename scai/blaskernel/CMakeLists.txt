###
 # @file scai/blaskernel/CMakeLists.txt
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
 # @brief scai/blaskernel/CMakeLists.txt
 # @author Eric Schricker
 # @date 11.11.2015
###

### Set PROJECT NAME / supported LANGUAGES
project ( scai_blaskernel )

### includes ###

## important Variables for project
set ( SCAI_BLASKERNEL_VERSION 1.0.0 )
include ( SetProjectVariables)
set ( SCAI_BLASKERNEL_EXTERNAL_DEPS   OpenMP CUDA MIC SCAI_BLAS                           ) 
include ( Dependencies/minimalSupportedVersions )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
include ( SetBuildFlags )
include ( VersionCheck/Compiler )
include ( Compiler/SetCPPFlags )

## Boost for Test
if ( BUILD_TEST )
    include( Package/BoostTest )
endif ( BUILD_TEST )

### include [optional] packages, use SCAI wrappers of CMake modules

include( Package/OpenMP )
include( Package/CUDA )
include( Package/MIC )
include( Package/SCAI_BLAS )

message( STATUS "CUDA: found = ${CUDA_FOUND}, inc = ${SCAI_CUDA_INCLUDE_DIR}, lib = ${SCAI_CUDA_LIBRARIES}" )
message( STATUS "SCAI_BLAS: found = ${SCAI_BLAS_FOUND}, name = ${SCAI_BLAS_NAME}, inc = ${SCAI_SCAI_BLAS_INCLUDE_DIR}, lib = ${SCAI_SCAI_BLAS_LIBRARIES}" )

### finalize build flags
include ( ConcludeFlags )

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

### Project System-Include ( also configured files )

set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )

include_directories ( ${SCAI_SCAI_BLAS_INCLUDE_DIR} )
include_directories ( ${SCAI_CUDA_INCLUDE_DIR} )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...
scai_project (

    CLASSES                  # .cpp, .hpp


    HEADERS                  # .hpp only

        BLASKernelTrait
        cblas
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

add_subdirectory ( doc )
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
add_library ( ${PROJECT_NAME} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )
set_target_properties ( ${PROJECT_NAME} PROPERTIES VERSION ${${UPPER_PROJECT_NAME}_VERSION} )

## add link libraries

target_link_libraries( ${PROJECT_NAME} scai_common scai_logging scai_tracing scai_tasking 
                                       scai_hmemo scai_kregistry )

target_link_libraries( ${PROJECT_NAME} ${SCAI_CUDA_LIBRARIES} )
target_link_libraries( ${PROJECT_NAME} ${SCAI_SCAI_BLAS_LIBRARIES} )

if ( USE_MIC )
    set ( MKLROOT $ENV{MKLROOT} )
    set ( CMAKE_SHARED_LINKER_FLAGS ${CMAKE_SHARED_LINKER_FLAGS} "-offload-option,mic,compiler,\"-L${MKLROOT}/lib/mic -Wl,--start-group -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -Wl,--end-group\"" )
    #ToDo: set_target_properties not working, any better way?
endif ( USE_MIC )

### install ###

## install lib
install ( TARGETS ${PROJECT_NAME} DESTINATION lib )

## install headers
relative_install( FILES ${CXX_HEADERS} DESTINATION ${${UPPER_PROJECT_NAME}_INCLUDE_DIR} )

## install library.hpp
install ( FILES ../${PROJECT_SURNAME}.hpp DESTINATION include/scai )

### add custum target distClean
include ( CustomCommands/DistClean )

