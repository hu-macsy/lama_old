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
include ( LAMAdefines ) # LAMA_: LOG_LEVEL, ASSERT_LEVEL, TRACE_LEVEL
include ( CompilerVersion )

include ( Package/OpenMP )
include ( Package/CUDA )
include ( Package/MIC )
include ( Package/Boost )
include ( Package/CodeCoverage ) # LAMA_USE_CODE_COVERAGE --> COVERAGE_FLAGS

## Boost 
if    ( Boost_INCLUDE_DIR )
    include_directories ( ${Boost_INCLUDE_DIR} )
endif ( Boost_INCLUDE_DIR )

## general compiler settings

include ( SetBuildFlags )
include ( SetCPPFlags )

## pass CACHE Variables to all subprojects

#get_cmake_property ( CACHE_VARS CACHE_VARIABLES )
#message ( STATUS "CACHE_VARS ${CACHE_VARS}" )

foreach    ( CACHE_VAR ${CACHE_VARS} )
    #get_property ( CACHE_VAR_HELPSTRING CACHE ${CACHE_VAR} PROPERTY HELPSTRING )
    #if    ( CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line." )
        
        get_property ( CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE )
        #message ( STATUS "CACHE_VAR ${CACHE_VAR}" )
        if     ( CACHE_VAR_TYPE STREQUAL "INTERNAL" )
            # skip Variable
        elseif ( CACHE_VAR_TYPE STREQUAL "UNINITIALIZED" )
            set ( CACHE_VAR_TYPE )
            list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
        else   ( )
            set( CACHE_VAR_TYPE :${CACHE_VAR_TYPE} )
            list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
        endif  ( )
        #set(CMAKE_ARGS "${CMAKE_ARGS} -D${CACHE_VAR}${CACHE_VAR_TYPE}=\"${${CACHE_VAR}}\"")
        #list ( APPEND CMAKE_ARGS "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
    #endif ( CACHE_VAR_HELPSTRING STREQUAL "No help, variable specified on the command line." )
endforeach ( CACHE_VAR ${CACHE_VARS} )

#message ( STATUS "Passed CMAKE_ARGS ${CMAKE_ARGS}" )

include ( ExternalProject )

ExternalProject_Add ( common 
            PREFIX common
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/common
            #CMAKE_ARGS "${CMAKE_ARGS}"
	    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
            )
                       
ExternalProject_Add ( logging 
            DEPENDS common
            PREFIX logging
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/logging
            #CMAKE_ARGS "${CMAKE_ARGS}"
	    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
            )
            
ExternalProject_Add ( tracing 
            DEPENDS common logging
            PREFIX tracing
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/tracing
            #CMAKE_ARGS "${CMAKE_ARGS}"
	    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
            )

ExternalProject_Add ( tasking 
            DEPENDS common logging tracing
            PREFIX tasking
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/tasking
            #CMAKE_ARGS "${CMAKE_ARGS}"
	    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
            )

ExternalProject_Add ( memory
            DEPENDS common logging tasking
            PREFIX memory
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/memory
            #CMAKE_ARGS "${CMAKE_ARGS}"
	    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
            )

ExternalProject_Add ( lama
            DEPENDS common logging memory tasking
            PREFIX lama
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/lama
            #CMAKE_ARGS "${CMAKE_ARGS}"
	    CMAKE_ARGS "-DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}"
            )
