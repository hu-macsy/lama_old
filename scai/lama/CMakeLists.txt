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

### Set PROJECT NAME / supported LANGUAGES
project ( LAMA CXX )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../CMake/Modules" )

include ( VersionDefinition )

### Include Compiler Flags
if    ( NOT SCAI_COMPLETE_BUILD )
	include ( SetBuildFlags )
	include ( SetCPPFlags )
endif ( NOT SCAI_COMPLETE_BUILD )

### include all search packages
include ( InternalDependencies )
include ( ExternalDependencies )

foreach    ( module ${SCAI_LAMA_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_LAMA_EXTERNAL_DEPS} )

foreach    ( PACKAGE_TO_FIND ${SCAI_LAMA_INTERNAL_DEPS} )
    find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND $${SCAI_LAMA_INTERNAL_DEPS} )

## finalize build flags
include ( ConcludeFlags )

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

### Project System-Include
include_directories ( ../.. )

foreach    ( module ${SCAI_LAMA_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    include_directories( ${${upper_module}_INCLUDE_DIR} )
endforeach ( module ${SCAI_LAMA_INTERNAL_DEPS} )

#foreach    ( module ${SCAI_LAMA_EXTERNAL_DEPS} )
#	string ( TOUPPER ${module} upper_module )
#    include_directories( ${${upper_module}_INCLUDE_DIR} )
#    message ( "add include dir for module ${module} to ${${upper_module}_INCLUDE_DIRS}" )
#endforeach ( module ${SCAI_LAMA_EXTERNAL_DEPS} )

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

include_directories ( ${SCAI_BLAS_INCLUDE_DIR} )

include ( SCAIAssert )

### build lama library

include ( Functions/lamaSourceSolution )

### set this dir as lama source dir
lama_set_source_dir ()

### Add classes and header files
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

if    ( NOT SCAI_BLAS_NAME STREQUAL "MKL" AND NOT INTERNALBLAS_FOUND )
    add_subdirectory( cblas )
endif ( NOT SCAI_BLAS_NAME STREQUAL "MKL" AND NOT INTERNALBLAS_FOUND )

## CUDA
if ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
    add_subdirectory ( cuda )
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )

## MIC
if ( USE_MIC )
    add_subdirectory ( mic )
endif ( USE_MIC )

add_subdirectory ( distribution )
add_subdirectory ( expression )
add_subdirectory ( io )
add_subdirectory ( macros )
add_subdirectory ( matrix )
add_subdirectory ( matutils )

## MPI
if ( MPI_FOUND AND USE_MPI )
    add_subdirectory ( mpi )
endif ( MPI_FOUND AND USE_MPI )

## GPI
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
	link_directories ( ${LAMAPP_LINK_DIRECTORIES} )
endif ( WIN32 )

### add library ###

if    ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )
	cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
	set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
endif ( CUDA_FOUND AND USE_CUDA AND CUDA_HAVE_GPU )

add_library ( ama ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} ) #${HEADERS} )

## add link libraries
foreach    ( module ${SCAI_LAMA_INTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( ama ${${upper_module}_LIBRARY} )
endforeach ( module ${SCAI_LAMA_INTERNAL_DEPS} )

foreach    ( module ${SCAI_LAMA_EXTERNAL_DEPS} )
	string ( TOUPPER ${module} upper_module )
    target_link_libraries ( ama ${SCAI_${upper_module}_LIBRARIES} )
endforeach ( module ${SCAI_LAMA_EXTERNAL_DEPS} )

# external libraries

if    ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )
    target_link_libraries ( ama ${GPI2_LIBRARIES} ${IBVERBS_LIBRARIES} )
endif ( GPI2_FOUND AND IBVERBS_FOUND AND USE_GPI2 )

if ( SCALAPACK_FOUND AND MPI_FOUND AND USE_MPI )
    target_link_libraries ( ama ${SCAI_PBLAS_LIBRARIES} )
    add_definitions( -DSCAI_MKL_SCALAPACK ) ## TODO: where used?!
endif ( SCALAPACK_FOUND AND MPI_FOUND AND USE_MPI )

### install
install ( TARGETS ama DESTINATION lib )

## install library.hpp
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/lama.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" DESTINATION include/scai )

# generate file make.inc + install it
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/examples/make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/examples/make.inc" DESTINATION share/examples )

# copy examples in build directory + install it
file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION ${CMAKE_CURRENT_BINARY_DIR} )
install ( DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/examples DESTINATION share )

## add custum target distClean, but only if this project is not included
if    ( NOT SCAI_COMPLETE_BUILD )
	include ( DistClean )
endif ( NOT SCAI_COMPLETE_BUILD )

#include ( LAMASummary )
