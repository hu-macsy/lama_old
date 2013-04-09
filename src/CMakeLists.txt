###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2011
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
 # $Id$
###

cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES and Module Path ###

project ( LAMA CXX C )
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

### Print LAMA Insallation header
message ( STATUS "" )
message ( STATUS "LAMA Installation:" )
message ( STATUS "==================" )
message ( STATUS "" )

### Variables / CMAKE Functions ###
include ( Variables )
include ( Functions )


if ( LAMA_CMAKE_VERBOSE )
    message ( STATUS "break 1" )
endif()

### PACKAGES ###
include ( PackageMPI )
include ( PackageBoost )
include ( PackageCodeCoverage )

# Find required packages
foreach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} QUIET REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
foreach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} QUIET )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )


### Include Compiler Flags, after all packages are found ###
include ( CompilerFlags )

## OPENMP
if ( NOT LAMA_OMP_SCHEDULE )
    set ( LAMA_OMP_SCHEDULE "static" )
endif ( NOT LAMA_OMP_SCHEDULE )

# TODO: SUMMARY
#message ( STATUS "OpenMP schedule set to " ${LAMA_OMP_SCHEDULE} )
add_definitions ( -DLAMA_OMP_SCHEDULE=${LAMA_OMP_SCHEDULE} )

if ( OPENMP_FOUND )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif( OPENMP_FOUND )




# ALLOW to switch off CUDA explicitly
setAndCheckCache ( "CUDA" )
setAndCheckCache ( "OPENCL" )

# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD )


## Project System-Include
include_directories ( . )

if ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIR} )
endif ( Boost_INCLUDE_DIR )


### add subdirectories ###
add_subdirectory ( lama )
add_subdirectory ( logging )

if ( LAMA_BUILD_TEST )
    add_subdirectory ( test )
endif ( LAMA_BUILD_TEST )

include ( PackageDoxygen )


### install ###
configure_file( "${CMAKE_SOURCE_DIR}/lama/lama.hpp.in" "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
# TODO: SUMMARY
#message ( STATUS "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )
install ( FILES "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" DESTINATION include )


### DISTCLEAN TARGET ###
add_custom_target ( distclean )
file ( GLOB_RECURSE BUILD_GLOB_RES ${CMAKE_BINARY_DIR}/* )
add_custom_command (
    TARGET distclean
    DEPENDS clean
    COMMAND ${CMAKE_COMMAND} -E remove ${BUILD_GLOB_RES}
)












### Summary ###

message( STATUS "" )
message( STATUS "" )
message( STATUS "Summary of LAMA Installation:" )
message( STATUS "=============================" )
message( STATUS "" )

message ( STATUS "Build Type    : ${CMAKE_BUILD_TYPE}" )

if ( LAMA_BLAS_LIBRARIES )
    message ( STATUS "BLAS libraries: ${LAMA_BLAS_LIBRARIES}" )
else ( LAMA_BLAS_LIBRARIES )
    message ( STATUS "BLAS not found" )
endif ( LAMA_BLAS_LIBRARIES )

if ( Boost_INCLUDE_DIR )
    message ( STATUS "Boost version ${Boost_VERSION} found at ${Boost_PATH}" )
    foreach ( lib ${Boost_COMPONENTS} )
        string ( TOUPPER ${lib} libname )
        set ( libname "Boost_${libname}_LIBRARY" ) 
        message ( STATUS "    ${libname} = ${${libname}}" )
    endforeach ( lib ${Boost_COMPONENTS} )
else ( Boost_INCLUDE_DIR )
    message ( STATUS "Boost not found" )
endif ( Boost_INCLUDE_DIR )
 
if( OPENMP_FOUND )
    message ( STATUS "OpenMP supported, uses compiler flags ${OpenMP_CXX_FLAGS}" )
else ( OPENMP_FOUND )
    message ( STATUS "OpenMP not found" )
endif ( OPENMP_FOUND )

if ( CUDA_FOUND )
    message ( STATUS "CUDA found at ${CUDA_INCLUDE_DIRS}" )
else ( CUDA_FOUND )
    message ( STATUS "CUDA not found" )
endif ( CUDA_FOUND )

if ( MPI_FOUND )
    message ( STATUS "MPI found at ${MPI_INCLUDE_PATH}" )
else ( MPI_FOUND )
    message ( STATUS "MPI not found" )
endif ( MPI_FOUND )

if ( OSHMEM_FOUND )
    message ( STATUS "OSHMEM found at ${OSHMEM_LIBRARIES}" )
else ( OSHMEM_FOUND )
    message ( STATUS "OSHMEM not found" )
endif ( OSHMEM_FOUND )

if ( GPI_FOUND )
    message ( STATUS "GPI found at ${GPI_LIBRARIES}" )
else ( GPI_FOUND )
    message ( STATUS "GPI not found" )
endif ( GPI_FOUND )

if ( GASNET_FOUND )
    message ( STATUS "GASNET found at ${GASNET_LIBRARIES}" )
else ( GASNET_FOUND )
    message ( STATUS "GASNET not found" )
endif ( GASNET_FOUND )

message ( STATUS "" )
message ( STATUS "LOG Level    : ${LAMA_LOG_LEVEL}" )
message ( STATUS "ASSERT Level : ${LAMA_ASSERT_LEVEL}" )
message ( STATUS "TRACE Level  : ${LAMA_TRACE_LEVEL}" )
