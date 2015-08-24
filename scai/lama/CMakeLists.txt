###
 # @file scai/lama/CMakeLists.txt
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
 # @brief Main specification file for CMake
 # @author Thomas Brandes, Jiri Kraus
 # @date 29.03.2012
 # @since 1.0.0
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

### Print LAMA Installation header
message ( STATUS "" )
message ( STATUS "LAMA Configuration:" )
message ( STATUS "===================" )
message ( STATUS "" )

### Set PROJECT NAME / supported LANGUAGES
project ( LAMA CXX )

### set Module Path

set ( LAMA_VERSION_NAME "Alpamayo" )
set ( LAMA_COPYRIGHT_YEAR 2015 )
set ( LAMA_VERSION_MAJOR 1 )
set ( LAMA_VERSION_MINOR 0 )
set ( LAMA_VERSION_PATCH 1 )
set ( LAMA_VERSION_COUNT 3 )
set ( LAMA_VERSION ${LAMA_VERSION_MAJOR}.${LAMA_VERSION_MINOR}.${LAMA_VERSION_PATCH} )

### Variables / CMAKE Functions ###
if    ( NOT SCAI_COMPLETE_BUILD )

	include ( SetBuildFlags )
	include ( SetCPPFlags )
	include ( CompilerVersion )
	
endif ( NOT SCAI_COMPLETE_BUILD )

message ( STATUS "CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS}" )

include ( Functions/allFunctions )

include ( SCAIAssert )

## Define all search packages
include ( InternalDependencies )
include ( ExternalDependencies )

# Find required packages

foreach    ( module ${SCAI_LAMA_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_LAMA_EXTERNAL_DEPS} )

# Find required packages
foreach    ( PACKAGE_TO_FIND ${SCAI_LAMA_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND $${SCAI_LAMA_INTERNAL_DEPS} )

## Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_LAMA_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_LAMA_INTERNAL_DEPS} )

if    ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIR} )
endif ( Boost_INCLUDE_DIR )

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

## MKL FOUND
if     ( SCAI_BLAS_NAME STREQUAL "MKL" )
    add_definitions( -DSCAI_MKL_BLAS )
elseif ( INTERNALBLAS_FOUND )
    add_definitions( -DSCAI_INTERNAL_BLAS )
else   ( )
    add_definitions( -DADD_ )
endif  ( )

add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )
add_definitions ( -D${SCAI_LOGGING_FLAG} )
add_definitions ( -D${SCAI_TRACING_FLAG} )

if ( WIN32 )
	add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

## build lama library

### set this dir as lama source dir ###
lama_set_source_dir ()

if    ( MPI_FOUND AND USE_MPI )
	include_directories ( ${MPI_INCLUDE_PATH} ) 
endif ( MPI_FOUND AND USE_MPI )

if    ( CUDA_FOUND AND USE_CUDA )
	include_directories ( ${CUDA_INCLUDE_DIRS} )
endif ( CUDA_FOUND AND USE_CUDA )

if ( METIS_FOUND )
    include_directories ( ${METIS_INCLUDE_DIR} ) 
endif ( METIS_FOUND )

if ( METIS_FOUND AND PARMETIS_FOUND )
    include_directories ( ${PARMETIS_INCLUDE_DIR} ) 
endif ( METIS_FOUND AND PARMETIS_FOUND )

if ( SCAI_BLAS_NAME STREQUAL "MKL" )
    include_directories ( ${MKL_INCLUDE_DIRS} )
elseif ( NOT SCAI_BLAS_NAME STREQUAL "MKL" AND NOT INTERNALBLAS_FOUND )
    include_directories( cblas/include )
#    add_subdirectory( cblas )
endif ( SCAI_BLAS_NAME STREQUAL "MKL" )

### Add classes and header files ###
set ( CLASSES
        CommunicationPlan
        Communicator
        Complex
        DenseVector
        Distributed
        LAMAArrayUtils
#        LAMAArrayView
        LAMAInterface
        LAMAInterfaceRegistry
        NoCommunicator
        Scalar
        StorageIO
        Vector
    )

set ( HEADERS
        LAMATypes
        TypeTraits
        BLASInterface
        UtilsInterface
        CRTPCommunicator
    )

lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )

### add subdirectories ###

if ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )

if ( USE_MIC )
    add_subdirectory ( mic )
endif ( USE_MIC )

add_subdirectory ( distribution )
add_subdirectory ( expression )
add_subdirectory ( io )
add_subdirectory ( macros )
add_subdirectory ( matrix )
add_subdirectory ( matutils )
if ( MPI_FOUND AND USE_MPI )
    add_subdirectory ( mpi )
endif ( MPI_FOUND AND USE_MPI )
if ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )
    add_subdirectory( gpi )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )
add_subdirectory ( norm )
add_subdirectory ( openmp )
add_subdirectory ( solver )
add_subdirectory ( storage )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

### set link directories ###
if ( WIN32 )
    set ( LAMAPP_LINK_DIRECTORIES ${LAMAPP_LINK_DIRECTORIES} ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

link_directories ( ${LAMAPP_LINK_DIRECTORIES} )

### add library ###

if    ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
	cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
	set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )

add_library ( ama ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} ) #${HEADERS} )

# internal libraries
foreach    ( module ${SCAI_LAMA_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( ama ${${upper_module}_LIBRARY} )
endforeach ( module ${SCAI_LAMA_INTERNAL_DEPS} )

# external libraries
target_link_libraries ( ama ${SCAI_BLAS_LIBRARIES} )

if    ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
	 target_link_libraries ( ama ${CUDA_cusparse_LIBRARY} ${CUDA_cublas_LIBRARY}
	 							 ${CUDA_CUDART_LIBRARY} ${CUDA_CUDA_LIBRARY} )          
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )

if    ( MPI_FOUND AND USE_MPI )
    target_link_libraries ( ama ${MPI_LIBRARIES} )
endif ( MPI_FOUND AND USE_MPI )

if    ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )
    target_link_libraries ( ama ${GPI2_LIBRARIES} ${IBVERBS_LIBRARIES} )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )

if ( SCALAPACK_FOUND AND MPI_FOUND AND USE_MPI )
    target_link_libraries ( ama ${SCAI_PBLAS_LIBRARIES} )
    add_definitions( -DSCAI_MKL_SCALAPACK ) ## TODO: where used?!
endif ( SCALAPACK_FOUND AND MPI_FOUND AND USE_MPI )

if ( METIS_FOUND )
    target_link_libraries( ama ${METIS_LIBRARY} )
endif ( METIS_FOUND ) 

if ( METIS_FOUND AND PARMETIS_FOUND )
    target_link_libraries( ama ${PARMETIS_LIBRARY} )
endif ( METIS_FOUND AND PARMETIS_FOUND )

# generate file make.inc + install it
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/examples/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" DESTINATION share/examples )

### add definitions ###
if ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### install ###

install ( TARGETS ama DESTINATION lib )

## end build lama library

### Examples ###
# copy examples in build directory + install it
file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION share )

## install library.hpp
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/lama.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" DESTINATION include/scai )

## add custum target distClean

## add custum target distClean, but only if this project is not included

if ( CMAKE_SOURCE_DIR EQUAL CMAKE_CURRENT_SOURCE_DIR )
    add_custom_target ( distclean )
endif ()

#include ( LAMASummary )
