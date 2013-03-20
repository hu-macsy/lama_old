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

##############################################################################
#
##############################################################################

cmake_minimum_required ( VERSION 2.8 )

### PROJECT NAME and supported LANGUAGES ###

project ( LAMA CXX C )

#### needed to find BLAS Libraries ###

if ( NOT WIN32 )
    enable_language ( Fortran )
endif ( NOT WIN32)

if ( CMAKE_C_COMPILER_ID MATCHES Intel )
    find_program ( XIAR xiar )
    if ( XIAR )
        set ( CMAKE_AR "${XIAR}" )
    endif ( XIAR )
    mark_as_advanced ( XIAR )
endif ( CMAKE_C_COMPILER_ID MATCHES Intel )

###  VARIABLES ###

# Set important variables
set ( CMAKE_SYSTEM_LIBRARY_PATH )
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )
set ( LAMA_ROOT_DIR "${CMAKE_SOURCE_DIR}/.." )

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

get_property ( FIND_LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS )
message ( STATUS "FindLib64: " ${FIND_LIB64} )

### COMPILER ###

include ( CompilerFlags )

### PACKAGES ###

# Look for MPI first to allow LAMA_BLAS to take the correct blacs implementation
# based on the found mpi
if ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )
    if ( MPI_ROOT )
        set ( MPI_C_INCLUDE_PATH "${MPI_ROOT}/Inc" )
        set ( MPI_CXX_INCLUDE_PATH "${MPI_ROOT}/Inc" )
        set ( LAMA_MPI_LIB_DIR "${MPI_ROOT}/Lib" )
    else ( MPI_ROOT )
        set ( MPI_C_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 R2/Inc" )
        set ( MPI_CXX_INCLUDE_PATH "C:/Program Files/Microsoft HPC Pack 2008 R2/Inc" )
        set ( LAMA_MPI_LIB_DIR "C:/Program Files/Microsoft HPC Pack 2008 R2/Lib" )
    endif ( MPI_ROOT )
    if ( CMAKE_CL_64 )
        set ( LAMA_MPI_LIB_DIR "${LAMA_MPI_LIB_DIR}/amd64" )
    else ( CMAKE_CL_64 )
        set ( LAMA_MPI_LIB_DIR "${LAMA_MPI_LIB_DIR}/i386" )
    endif ( CMAKE_CL_64 )
    set ( MPI_C_LIBRARIES "${LAMA_MPI_LIB_DIR}/msmpi.lib" )
    set ( MPI_CXX_LIBRARIES "${LAMA_MPI_LIB_DIR}/msmpi.lib" )
else ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) )
  if ( MPI_ROOT )
      set ( MPI_COMPILER ${MPI_ROOT}/bin/mpicxx )
  endif ( MPI_ROOT )
endif ( WIN32 AND NOT ( MPI_C_INCLUDE_PATH OR MPI_CXX_INCLUDE_PATH OR MPI_C_LIBRARIES OR MPI_CXX_LIBRARIES ) ) 

##############################################################################
#  MPI Stuff
##############################################################################

find_package ( MPI )

### ALLOW to switch off MPI explicitly ###

if ( NOT LAMA_USE_MPI )
    set ( DEFAULT_USE_MPI ${MPI_FOUND} )
endif ( NOT LAMA_USE_MPI )

set ( LAMA_USE_MPI ${DEFAULT_USE_MPI} CACHE BOOL "Enable / Disable use of MPI" )

if ( NOT LAMA_USE_MPI )
    set ( MPI_FOUND FALSE )
endif ( NOT LAMA_USE_MPI )

message( STATUS "MPI: found = ${MPI_FOUND}, use = ${LAMA_USE_MPI}" )

if ( MPI_FOUND )
    include_directories ( ${MPI_INCLUDE_PATH} )
endif ( MPI_FOUND )

##############################################################################
#  Code coverage with gcov/lcov
##############################################################################

if ( NOT CODE_COVERAGE )
    set ( DEFAULT_CODE_COVERAGE FALSE )
endif ( NOT CODE_COVERAGE )

set ( CODE_COVERAGE ${DEFAULT_CODE_COVERAGE} CACHE BOOL "Enable / Disable Code Coverage" )

if ( CODE_COVERAGE )
	set ( COVERAGE_FLAGS "-fprofile-arcs -ftest-coverage" )
	set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COVERAGE_FLAGS}" )
endif ( CODE_COVERAGE ) 

##############################################################################
#  Other packages
##############################################################################

# Find required packages
set ( REQUIRED_PACKAGES_TO_FIND
        LAMA_BLAS
        #add required packages here
    ) 

foreach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} QUIET REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
set ( OPTIONAL_PACKAGES_TO_FIND
        OpenMP
        Doxygen
        Threads
        #DirectX
        #OSHMEM REQUIRES GASNET
        #GPI
        #add optional packages here
    )

#CUDA Only works with GCC on Linux
#TODO: This needs to be checked on windows
if ( CMAKE_COMPILER_IS_GNUCC )
    set ( OPTIONAL_PACKAGES_TO_FIND
          ${OPTIONAL_PACKAGES_TO_FIND}
          CUDA
	      #CUDA needs to be found before OpenCL because CUDA_INCLUDE_DIRS is used as a hint to find the OpenCL Headers
    	  #OpenCL
    )
endif ( CMAKE_COMPILER_IS_GNUCC )

if ( NOT DEFINED DISABLE_METIS )
  set ( OPTIONAL_PACKAGES_TO_FIND
          ParMetis REQUIRES Metis
          ${OPTIONAL_PACKAGES_TO_FIND} )
endif ( NOT DEFINED DISABLE_METIS )

foreach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} QUIET )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )

### BOOST Specific ###

# set ( Boost_USE_MULTITHREADED OFF )

if ( WIN32 )
    message ( STATUS "Setting special Boost options on Windows" )
    #set ( Boost_USE_STATIC_LIBS ON )
    set ( Boost_USE_MULTITHREADED ON )
    set ( Boost_USE_STATIC_RUNTIME OFF )
    #add_definitions ( -DBOOST_ALL_NO_LIB )
endif ( WIN32 )

# Finds packages with custom search options 

set ( Boost_COMPONENTS thread program_options unit_test_framework regex )

# FindBoost Debug options comment
if ( LAMA_DEBUG_CMAKE )
    set ( Boost_DEBUG TRUE )
    set ( Boost_DETAILED_FAILURE_MSG TRUE )
endif( LAMA_DEBUG_CMAKE )

# Find Boost 

find_package ( Boost COMPONENTS ${Boost_COMPONENTS} QUIET )

# Note: we use Boost_INCLUDE_DIR, Boost_<lib>_FOUND, Boost_<lib>_LIBRARY, but
#       not Boost_FOUND, as it is false if some optional libraries are missing

# Boost: include directory is mandatory ( LAMA uses shared pointer, function )

if ( Boost_INCLUDE_DIR )
    message ( STATUS "Boost_INCLUDE_DIR = ${Boost_INCLUDE_DIR}" )
    get_filename_component ( Boost_PATH ${Boost_INCLUDE_DIR} PATH )
    message ( STATUS "Boost_PATH = ${Boost_PATH}" )
    # Boost_PATH should be same as BOOST_ROOT
else ( Boost_INCLUDE_DIR )
    message ( FATAL_ERROR "Boost (include directory) not found: give hint by environment variable BOOST_ROOT" ) 
endif ( Boost_INCLUDE_DIR )

# check status of each Boost component

foreach ( lib ${Boost_COMPONENTS} )
   string ( TOUPPER ${lib} libname )
   set ( libname "Boost_${libname}_LIBRARY" )
   # libname: variable that contains the library for the boost component
   message ( STATUS "${libname} = ${${libname}}" )
   if ( ${libname} )
       # library found, make sure it belongs to same version of Boost
       if ( "${${libname}}" MATCHES "${Boost_PATH}*" )
       	   #
       else ( "${${libname}}" MATCHES "${Boost_PATH}*" )
           message ( FATAL_ERROR "${${libname}} illegal, not in ${Boost_PATH}" )
       endif ( "${${libname}}" MATCHES "${Boost_PATH}*" )
   endif ( ${libname} )
endforeach ( lib ${Boost_COMPONENTS} )

if ( NOT Boost_THREAD_LIBRARY )
    message ( FATAL_ERROR "Boost thread library not found" )
endif ( NOT Boost_THREAD_LIBRARY )

## MAKEFILES

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF ) 

## BUILDTYPE

# Choose Default CMAKE_BUILD_TYPE
if ( NOT CMAKE_BUILD_TYPE )
  # Can be: (RelWithDebInfo)
  set ( CMAKE_BUILD_TYPE Debug CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

message ( STATUS "Build type is set to " ${CMAKE_BUILD_TYPE} )

## LOGGING Level
#
#  Debug   : use -DLOG_LEVEL_DEBUG
#  Release : use -DLOG_LEVEL_INFO
#  
#  For serious problems: -DLOG_LEVEL_TRACE
#  For benchmarks:       -DLOG_LEVEL_OFF (or -DLOG_LEVEL_FATAL, -DLOG_LEVEL_ERROR)

if ( NOT LAMA_LOG_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_LOG_LEVEL "INFO" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    else ()
        set ( DEFAULT_LOG_LEVEL "TRACE" )
    endif ()
endif ( NOT LAMA_LOG_LEVEL )

set ( LAMA_LOG_LEVEL ${DEFAULT_LOG_LEVEL} CACHE STRING
      "Choose level of compile time logging: TRACE, DEBUG, INFO, WARN, ERROR, OFF" )

add_definitions ( -DLAMA_LOG_LEVEL_${LAMA_LOG_LEVEL} )

## ASSERT Level
#
#  Debug   : use -DASSERT_LEVEL_DEBUG
#  Release : use -DASSERT_LEVEL_ERROR
#  
#  For benchmarks:       -DASSERT_LEVEL_OFF

if ( NOT LAMA_ASSERT_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_ASSERT_LEVEL "ERROR" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    else ()
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    endif ()
endif ( NOT LAMA_ASSERT_LEVEL )

set ( LAMA_ASSERT_LEVEL ${DEFAULT_ASSERT_LEVEL} CACHE STRING
      "Choose level of ASSERT: DEBUG, ERROR, OFF" )

add_definitions ( -DLAMA_ASSERT_LEVEL_${LAMA_ASSERT_LEVEL} )

## LAMA TRACE LEVEL
#
# If TRACE is set to OFF all LAMA_REGION macros in the code are
# completely ignored. If TRACE is set to VT, regions will be traced
# (entry, exit event) for VampirTrace.

if ( NOT LAMA_TRACE_LEVEL )
    set ( DEFAULT_TRACE_LEVEL "OFF" )
endif ( NOT LAMA_TRACE_LEVEL )

set ( LAMA_TRACE_LEVEL ${DEFAULT_TRACE_LEVEL} CACHE STRING
     "Choose level of TRACE: VT (for VampirTrace), TIME(region timing), SIMPLE(simple timing) or OFF (default)" )

add_definitions( -DLAMA_TRACE_LEVEL_${LAMA_TRACE_LEVEL} )

## OPENMP
if ( NOT LAMA_OMP_SCHEDULE )
    set ( LAMA_OMP_SCHEDULE "static" )
endif ( NOT LAMA_OMP_SCHEDULE )

message ( STATUS "OpenMP schedule set to " ${LAMA_OMP_SCHEDULE} )
add_definitions ( -DLAMA_OMP_SCHEDULE=${LAMA_OMP_SCHEDULE} )

if ( OPENMP_FOUND )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif( OPENMP_FOUND )

## ALLOW to switch off CUDA explicitly

if ( NOT LAMA_USE_CUDA )
    set ( DEFAULT_USE_CUDA ${CUDA_FOUND} )
endif ( NOT LAMA_USE_CUDA )

set ( LAMA_USE_CUDA ${DEFAULT_USE_CUDA} CACHE BOOL "Enable / Disable use of CUDA" )

if ( NOT LAMA_USE_CUDA )
    set ( CUDA_FOUND FALSE )
endif ( NOT LAMA_USE_CUDA )

## ALLOW to switch off OPENCL explicitly

if ( NOT LAMA_USE_OPENCL )
    set ( DEFAULT_USE_OPENCL ${OPENCL_FOUND} )
endif ( NOT LAMA_USE_OPENCL )

set ( LAMA_USE_OPENCL ${DEFAULT_USE_OPENCL} CACHE BOOL "Enable / Disable use of OPENCL" )

if ( NOT LAMA_USE_OPENCL )
    set ( OPENCL_FOUND FALSE )
endif ( NOT LAMA_USE_OPENCL )

## CUDA SPECIFIC

# LAMA irrelevant entries will be marked as advanced

mark_as_advanced ( CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD )

## Project System-Include

include_directories ( . )

if ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIR} )
endif ( Boost_INCLUDE_DIR )

if ( OPENCL_FOUND )
    add_definitions ( -DLAMA_BUILD_OPENCL )
    include_directories ( ${OPENCL_INCLUDE_DIRS} )
endif ( OPENCL_FOUND )

## Library System-Includes

if ( OSHMEM_FOUND )
    message ( STATUS "OpenShMem found, include dir = ${OSHMEM_INCLUDE_DIR}" )
else ( OSHMEM_FOUND )
    message ( STATUS "OpenShMem not found" )
endif ( OSHMEM_FOUND )

if ( GPI_FOUND )
    message ( STATUS "GPI found, include dir = ${GPI_INCLUDE_DIR}" )
else ( GPI_FOUND )
    message ( STATUS "GPI not found" )
endif ( GPI_FOUND )


### add subdirectories ###

add_subdirectory ( lama )
add_subdirectory ( logging )

if ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )
    add_subdirectory ( test )
else ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )
    if ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND )
       message ( WARNING "Not building tests because Boost unit test framework is missing." )
    endif ( NOT Boost_UNIT_TEST_FRAMEWORK_FOUND )
    if ( NOT Boost_REGEX_FOUND )
       message ( WARNING "Not building tests because Boost regex is missing." )
    endif ( NOT Boost_REGEX_FOUND )
    if ( ${Boost_VERSION} LESS "104100" )
       message ( WARNING "Not building tests because Boost is to old: ${Boost_VERSION}." )
    endif ( ${Boost_VERSION} LESS "104100" )
endif ( ${Boost_VERSION} GREATER "104099" AND Boost_UNIT_TEST_FRAMEWORK_FOUND AND Boost_REGEX_FOUND )

### DOXYGEN DOCUMENTATION ###

if ( DOXYGEN_FOUND )

    set ( LAMA_DOC_DIR "${LAMA_SOURCE_DIR}/../doc" )
    set ( DOXYGEN_BUILD_ROOT "${LAMA_DOC_DIR}/doxygen" )
    filE ( MAKE_DIRECTORY ${DOXYGEN_BUILD_ROOT} )

   # The initial rm command gets rid of everything previously built by this
   # custom command.

   add_custom_command (
      OUTPUT ${DOXYGEN_BUILD_ROOT}/html/index.html
      COMMAND rm -rf ${DOXYGEN_BUILD_ROOT}
      COMMAND mkdir ${DOXYGEN_BUILD_ROOT}
      COMMAND ${DOXYGEN_EXECUTABLE} ${LAMA_DOC_DIR}/LAMA.Doxyfile
      DEPENDS ${LAMA_DOC_DIR}/LAMA.Doxyfile
      WORKING_DIRECTORY ${LAMA_DOC_DIR}
   )

   add_custom_target (
      doc
      DEPENDS
      ${DOXYGEN_BUILD_ROOT}/html/index.html
   )

   # Install the documentation generated at "make" time.

   # install ( DIRECTORY ${DOXYGEN_BUILD_ROOT}/ DESTINATION ${DOXYGEN_BUILD_ROOT}/html )

else ( DOXYGEN_FOUND )
    message ( WARNING "Not building system documentation because Doxygen not found." )
endif ( DOXYGEN_FOUND )

file ( MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/include")

### boost version check ###

set ( BOOST_VERSION_ERROR_STRING "Included Boost version does not match linked Boost version ${Boost_VERSION}" )

set ( LAMA_HPP_CONTENS " 
#include <boost/version.hpp>
#define LAMA_LINKED_BOOST_VERSION ${Boost_VERSION}

#if LAMA_LINKED_BOOST_VERSION != BOOST_VERSION
    #error \"${BOOST_VERSION_ERROR_STRING}\"
#endif //LAMA_LINKED_BOOST_VERSION != BOOST_VERSION

#define LAMA_REVISION ${Project_WC_REVISION}

#define LAMA_ASSERT_LEVEL_${LAMA_ASSERT_LEVEL}
#define LAMA_LOG_LEVEL_${LAMA_LOG_LEVEL}
#define LAMA_TRACE_LEVEL_${LAMA_TRACE_LEVEL}
#define LAMA_OMP_SCHEDULE ${LAMA_OMP_SCHEDULE}

extern \"C\" inline int lama_getRevision()
{
    return LAMA_REVISION\;
}

")


### install ###

file ( WRITE "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" ${LAMA_HPP_CONTENS} )

message ( STATUS "${CMAKE_CURRENT_BINARY_DIR}/include/lama.hpp" )

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

if ( OPENCL_FOUND )
    message ( STATUS "OpenCL found at ${OPENCL_INCLUDE_DIRS}" )
else ( OPENCL_FOUND )
    message ( STATUS "OpenCL not found" )
endif ( OPENCL_FOUND )

if ( MPI_FOUND )
    message ( STATUS "MPI found at ${MPI_INCLUDE_PATH}" )
else ( MPI_FOUND )
    message ( STATUS "MPI not found" )
endif ( MPI_FOUND )

if ( METIS_FOUND )
    message ( STATUS "Metis found at ${METIS_INCLUDE_DIR}" )
else( METIS_FOUND )
    message ( STATUS "Metis not found" )
endif( METIS_FOUND )

if ( PARMETIS_FOUND )
    message ( STATUS "ParMetis found at ${PARMETIS_INCLUDE_DIR}" )
else( PARMETIS_FOUND )
    message ( STATUS "ParMetis not found" )
endif( PARMETIS_FOUND )

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
