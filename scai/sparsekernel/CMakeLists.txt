###
 # @file scai/sparsekernel/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2016
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

### Set required CMake version
cmake_minimum_required ( VERSION 2.8.8 )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### set cmake quiet for untyped messages
if    ( NOT SCAI_CMAKE_VERBOSE )
    include ( Functions/setCMakeQuiet )
endif ( NOT SCAI_CMAKE_VERBOSE )

### Set PROJECT NAME / supported LANGUAGES
project ( scai_sparsekernel CXX )

### includes ###

## important Variables for project
include ( VersionDefinition )
include ( SetProjectVariables)
include ( Dependencies/internal )
include ( Dependencies/external )
include ( Dependencies/minimalSupportedVersions )

## needed functions
include ( Functions/scaiFunctions )
include ( Functions/scaiProject )
include ( Functions/relativeInstall )

### Include compile dependent things
if    ( NOT SCAI_COMPLETE_BUILD )
    include ( SetBuildFlags )
    include ( VersionCheck/Compiler )
    include ( Compiler/SetCPPFlags )
    include ( Package/Sphinx )
else  ( NOT SCAI_COMPLETE_BUILD )
    include ( Package/doc )
    include ( Functions/CheckValueAtCompleteBuild )
    CheckValueAtCompleteBuild ( ${PROJECT_SURNAME} )
endif ( NOT SCAI_COMPLETE_BUILD )

## Boost for Test
if ( BUILD_TEST )
    include( Package/BoostTest )
endif ( BUILD_TEST )

### include all search packages
addInternalAndExternalPackages()

### finalize build flags
include ( ConcludeFlags )

### add definitions
add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if ( SCAI_OMP_SCHEDULE_FLAG )
    add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
endif ( SCAI_OMP_SCHEDULE_FLAG )

add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -D${SCAI_TRACING_FLAG} )

if    ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

## Project System-Include
include_directories ( ../.. ${CMAKE_BINARY_DIR}/include )
addInternalAndExternalIncludes()

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
addInternalAndExternalLinkLibraries()

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

### summary ###
if    ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
    include ( Summaries/${PROJECT_SURNAME} )
endif ( DEFINED SCAI_CMAKE_VERBOSE OR NOT DEFINED SCAI_COMPLETE_BUILD )
