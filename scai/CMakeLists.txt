###
 # @file CMakeLists.txt
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
 # @author Lauretta Schubert
 # @date 30.07.2015
 # @since 2.0.0
###

cmake_minimum_required ( VERSION 2.8 )

project ( LAMA_ALL CXX )
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

include ( VersionDefinition )

## general Settings for all projects
include ( SetBuildFlags ) # VERBOSE, RPATH, LIB64
include ( SCAIAssert )
include ( SetCPPFlags )

include ( CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 )

include ( Functions/checkValue )
include ( Functions/setAndCheckCache )
include ( Functions/checkWhiteList )
include ( Functions/concatCacheArguments )

include ( ExternalDependencies )

foreach    ( module ${SCAI_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${SCAI_EXTERNAL_DEPS} )

include ( LoggingSettings )
include ( TracingSettings )

include ( ConcludeFlags )

## documentation

# user docu --> make userdoc, make userdoc_json, make install_userdoc
include ( Package/Sphinx )

#system docu --> make doc
include ( Package/Doxygen )

## summary of what is found
include ( Summaries/scai )

## pass CACHE Variables to all subprojects via white listing

# loads variables WHITELIST, MPI_WHITELIST, CUDA_WHITELIST, MIC_WHITELIST
include ( NeededCacheVariables )

set ( BOOST_ARGS "" )
set ( CUDA_ARGS "" )
set ( GPI2_ARGS "" )
set ( GRAPHPARTITIONING_ARGS "" )
set ( MIC_ARGS "" )
set ( MPI_ARGS "" )
set ( OPENMP_ARGS "" )
set ( SCAI_BLAS_ARGS "" )

get_cmake_property ( CACHE_VARS CACHE_VARIABLES )
foreach    ( CACHE_VAR ${CACHE_VARS} )

	#message ( STATUS "test VAR ${CACHE_VAR}" )
	check_whitelist ( CACHE_VAR "${WHITELIST}" "CMAKE_ARGS" )
	
	check_whitelist ( CACHE_VAR "${BOOST_WHITELIST}" "BOOST_ARGS" )
	check_whitelist ( CACHE_VAR "${CUDA_WHITELIST}" "CUDA_ARGS" )
	check_whitelist ( CACHE_VAR "${GPI2_WHITELIST}" "GPI2_ARGS" )
	check_whitelist ( CACHE_VAR "${GRAPHPARTIONING_WHITELIST}" "GRAPHPARTIONING_ARGS" )
	check_whitelist ( CACHE_VAR "${MIC_WHITELIST}" "MIC_ARGS" )    
	check_whitelist ( CACHE_VAR "${MPI_WHITELIST}" "MPI_ARGS" )
	check_whitelist ( CACHE_VAR "${OPENMP_WHITELIST}" "OPENMP_ARGS" )
	check_whitelist ( CACHE_VAR "${SCAI_BLAS_WHITELIST}" "SCAI_BLAS_ARGS" )
    
endforeach ( CACHE_VAR ${CACHE_VARS} )

#message ( STATUS "Passed CMAKE_ARGS ${CMAKE_ARGS}" )

#message ( STATUS "Passed BOOST_ARGS ${BOOST_ARGS}" )
#message ( STATUS "Passed CUDA_ARGS ${CUDA_ARGS}" )
#message ( STATUS "Passed GPI2_ARGS ${GPI_ARGS}" )
#message ( STATUS "Passed GRAPHPARTIONING_ARGS ${GRAPHPARTIONING_ARGS}" )
#message ( STATUS "Passed MIC_ARGS ${MIC_ARGS}" )
#message ( STATUS "Passed MPI_ARGS ${MPI_ARGS}" )
#message ( STATUS "Passed OPENMP_ARGS ${OPENMP_ARGS}" )
#message ( STATUS "Passed SCAI_BLAS_ARGS ${SCAI_BLAS_ARGS}" )

## define external projects

# loads SCAI_[LIB]_INTERNAL_DEPS variables 
include ( InternalDependencies )
include ( ExternalDependencies )

include ( ExternalProject )

#message ( STATUS "general ${CMAKE_ARGS}" )

## TODO: do in foreach loop
#foreache ( project scai_common scai_common scai_logging scai_tracing scai_tasking scai_hmemo scai_lama )

concatCacheArguments ( scai_common "CACHE_ARGS" )
#message ( STATUS "common CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_common 
			DEPENDS ${SCAI_COMMON_INTERNAL_DEPS}
            PREFIX common
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/common
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
	    	BINARY_DIR common
            )

concatCacheArguments ( scai_logging "CACHE_ARGS" )
#message ( STATUS "logging CACHE_ARGS ${CACHE_ARGS}" )
         
ExternalProject_Add ( scai_logging 
            DEPENDS ${SCAI_LOGGING_INTERNAL_DEPS}
            PREFIX logging
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/logging
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
	    	BINARY_DIR logging
            )

concatCacheArguments ( scai_tracing "CACHE_ARGS" )
#message ( STATUS "tracing CACHE_ARGS ${CACHE_ARGS}" )
           
ExternalProject_Add ( scai_tracing 
            DEPENDS ${SCAI_TRACING_INTERNAL_DEPS}
            PREFIX tracing
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/tracing
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
	    	BINARY_DIR tracing
            )

concatCacheArguments ( scai_tasking "CACHE_ARGS" )
#message ( STATUS "tasking CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_tasking 
            DEPENDS ${SCAI_TASKING_INTERNAL_DEPS}
            PREFIX tasking
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/tasking
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
	    	BINARY_DIR tasking
            )

concatCacheArguments ( scai_kernel "CACHE_ARGS" )
#message ( STATUS "kernel CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_kernel
            DEPENDS ${SCAI_HMEMO_INTERNAL_DEPS}
            PREFIX kernel
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/kernel
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
	    	BINARY_DIR kernel
            )

concatCacheArguments ( scai_lama "CACHE_ARGS" )
#message ( STATUS "lama CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_lama
            DEPENDS ${SCAI_LAMA_INTERNAL_DEPS}
            PREFIX lama
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/lama
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
	    	BINARY_DIR lama
            )

## add custom target docclean
include ( DocClean )

## add custom target distclean
include ( DistClean )
