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

## version info
set ( LAMA_VERSION_NAME "Alpamayo" )
set ( LAMA_COPYRIGHT_YEAR 2015 )
set ( LAMA_VERSION_MAJOR 1 )
set ( LAMA_VERSION_MINOR 0 )
set ( LAMA_VERSION_PATCH 1 )
set ( LAMA_VERSION_COUNT 3 )
set ( LAMA_VERSION ${LAMA_VERSION_MAJOR}.${LAMA_VERSION_MINOR}.${LAMA_VERSION_PATCH} )

## general Settings for all projects

include ( SetBuildFlags ) # VERBOSE, RPATH, LIB64, BUILD_SHARED
include ( SCAIAssert )
include ( CompilerVersion )

include ( Functions/checkValue )
include ( Functions/setAndCheckCache )
include ( Functions/checkWhiteList )

include ( ExternalDependencies )

foreach    ( module ${SCAI_EXTERNAL_DEPS} )
message ( STATUS "include package ${module}" )
    include( Package/${module} )
endforeach ( module ${SCAI_EXTERNAL_DEPS} )

include ( SetCPPFlags ) # do after external deps because openmp dependency

## documentation

# user docu --> make userdoc, make userdoc_json, make install_userdoc
include ( Package/Sphinx )

#system docu --> make doc
include ( Package/Doxygen )

## summary of what is found
include ( Summary )

## pass CACHE Variables to all subprojects via white listing

# loads variables WHITELIST, MPI_WHITELIST, CUDA_WHITELIST, MIC_WHITELIST
include ( NeededCacheVariables )

set ( MPI_ARGS "" )
set ( CUDA_ARGS "" )
set ( MIC_ARGS "" )

get_cmake_property ( CACHE_VARS CACHE_VARIABLES )
foreach    ( CACHE_VAR ${CACHE_VARS} )

	check_whitelist ( CACHE_VAR "${WHITELIST}" "${CMAKE_ARGS}" )
	check_whitelist ( CACHE_VAR "${MPI_WHITELIST}" "${MPI_ARGS}" )  
	check_whitelist ( CACHE_VAR "${CUDA_WHITELIST}" "${CUDA_ARGS}" )
	check_whitelist ( CACHE_VAR "${MIC_WHITELIST}" "${MIC_ARGS}" )    
    
endforeach ( CACHE_VAR ${CACHE_VARS} )

message ( STATUS "Passed CMAKE_ARGS ${CMAKE_ARGS}" )
#message ( STATUS "Passed MPI_ARGS ${MPI_ARGS}" )
#message ( STATUS "Passed CUDA_ARGS ${CUDA_ARGS}" )
#message ( STATUS "Passed MIC_ARGS ${MIC_ARGS}" )

## define external projects

# loads SCAI_[LIB]_INTERNAL_DEPS variables 
include ( InternalDependencies )

include ( ExternalProject )

ExternalProject_Add ( scai_common 
			DEPENDS ${SCAI_COMMON_INTERNAL_DEPS}
            PREFIX common
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/common
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS} ${CUDA_ARGS}"
	    	BINARY_DIR common
            )
                       
ExternalProject_Add ( scai_logging 
            DEPENDS ${SCAI_LOGGING_INTERNAL_DEPS}
            PREFIX logging
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/logging
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS}"
	    	BINARY_DIR logging
            )
            
ExternalProject_Add ( scai_tracing 
            DEPENDS ${SCAI_TRACING_INTERNAL_DEPS}
            PREFIX tracing
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/tracing
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS}"
	    	BINARY_DIR tracing
            )

ExternalProject_Add ( scai_tasking 
            DEPENDS ${SCAI_TASKING_INTERNAL_DEPS}
            PREFIX tasking
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/tasking
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS}"
	    	BINARY_DIR tasking
            )

ExternalProject_Add ( scai_hmemo
            DEPENDS ${SCAI_HMEMO_INTERNAL_DEPS}
            PREFIX hmemo
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/hmemo
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS} ${CUDA_ARGS}" #${MIC_ARGS}
	    	BINARY_DIR hmemo
            )

ExternalProject_Add ( lama
            DEPENDS ${SCAI_LAMA_INTERNAL_DEPS}
            PREFIX lama
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/lama
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS} ${CUDA_ARGS}" #${MIC_ARGS} ${MPI_ARGS}
	    	BINARY_DIR lama
            )
