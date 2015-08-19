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

project ( LAMA_ALL CXX C )
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

## general Settings for all projects

include ( Variables ) # VERBOSE, RPATH, LIB64, BUILD_SHARED
include ( Functions/checkValue )
include ( Functions/setAndCheckCache )
include ( SCAIdefines ) # SCAI_: LOG_LEVEL, ASSERT_LEVEL, TRACE_LEVEL
include ( CompilerVersion )

include ( Package/OpenMP )
include ( Package/CUDA )
include ( Package/MIC )
include ( Package/Boost )
#include ( Package/MPI )
include ( Package/CodeCoverage ) # LAMA_USE_CODE_COVERAGE --> COVERAGE_FLAGS

set ( REQUIRED_PACKAGES_TO_FIND Threads )
set ( OPTIONAL_PACKAGES_TO_FIND LAMA_BLAS ParMetis REQUIRES Metis )

# Find required packages
foreach    ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} REQUIRED )
endforeach ( PACKAGE_TO_FIND ${REQUIRED_PACKAGES_TO_FIND} )

# Find optional packages
foreach    ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )
    find_package ( ${PACKAGE_TO_FIND} ${LAMA_FIND_PACKAGE_FLAGS} )
endforeach ( PACKAGE_TO_FIND ${OPTIONAL_PACKAGES_TO_FIND} )

## documentation

# user docu --> make userdoc, make userdoc_json, make install_userdoc
find_package ( Sphinx )
include ( Package/Sphinx )

#system docu --> make doc
find_package ( Doxygen )
include ( Package/Doxygen )

## general compiler settings

include ( SetBuildFlags )
include ( SetCPPFlags )
include ( SetNVCCFlags )


## summary of what is found



## pass CACHE Variables to all subprojects via white listing

# loads variables WHITELIST, CUDA_WHITELIST, MPI_WHITELIST, MIC_WHITELIST
include ( NeededCacheVariables )

get_cmake_property ( CACHE_VARS CACHE_VARIABLES )
#message ( STATUS "CACHE_VARS ${CACHE_VARS}" )

foreach    ( CACHE_VAR ${CACHE_VARS} )
    foreach    ( WHITELIST_ITEM ${WHITELIST} )
        
        string( COMPARE EQUAL "${CACHE_VAR}" "${WHITELIST_ITEM}" result )
        if    ( ${result} )
        	#message ( STATUS "${CACHE_VAR} matches ${WHITELIST_ITEM}" )
        
	        get_property ( CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE )
	        #message ( STATUS "CACHE_VAR ${CACHE_VAR}" )
	        if     ( CACHE_VAR_TYPE STREQUAL "INTERNAL" )
	            # skip Variable
	        #elseif ( CACHE_VAR_TYPE STREQUAL "UNINITIALIZED" )
	        #    set ( CACHE_VAR_TYPE )
	        #    list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        else   ( )
	            set( CACHE_VAR_TYPE :${CACHE_VAR_TYPE} )
	            list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        endif  ( )
	        #set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
	        #list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
        
        endif ( ${result} )
        
    endforeach ( WHITELIST_ITEM ${WHITELIST} )
    
 foreach    ( WHITELIST_ITEM ${CUDA_WHITELIST} )
        
        string( COMPARE EQUAL "${CACHE_VAR}" "${WHITELIST_ITEM}" result )
        if    ( ${result} )
        	#message ( STATUS "${CACHE_VAR} matches ${WHITELIST_ITEM}" )
        
	        get_property ( CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE )
	        #message ( STATUS "CACHE_VAR ${CACHE_VAR}" )
	        if     ( CACHE_VAR_TYPE STREQUAL "INTERNAL" )
	            # skip Variable
	        #elseif ( CACHE_VAR_TYPE STREQUAL "UNINITIALIZED" )
	        #    set ( CACHE_VAR_TYPE )
	        #    list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        else   ( )
	            set( CACHE_VAR_TYPE :${CACHE_VAR_TYPE} )
	            list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        endif  ( )
	        #set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
	        #list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
        
        endif ( ${result} )
        
    endforeach ( WHITELIST_ITEM ${CUDA_WHITELIST} )
    
    ## do the same for MPI_WHITELIST, MIC_WHITELIST
    
endforeach ( CACHE_VAR ${CACHE_VARS} )

#message ( STATUS "Passed CMAKE_ARGS ${CMAKE_ARGS}" )

## define external projects

# loads SCAI_[LIB]_INTERNAL_DEPS variables 
include ( InternalDependencies )

include ( ExternalProject )

ExternalProject_Add ( scai_common 
			DEPENDS ${SCAI_COMMON_INTERNAL_DEPS}
            PREFIX common
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/common
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS}"
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
            CMAKE_CACHE_ARGS "${CMAKE_ARGS}"
	    	BINARY_DIR hmemo
            )

ExternalProject_Add ( lama
            DEPENDS ${SCAI_LAMA_INTERNAL_DEPS}
            PREFIX lama
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/lama
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS}"
	    	BINARY_DIR lama
            )
