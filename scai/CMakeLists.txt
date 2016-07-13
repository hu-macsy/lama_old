###
 # @file CMakeLists.txt
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
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 30.07.2015
###

### Set required CMake version
cmake_minimum_required ( VERSION 2.8.8 )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../CMake/Modules" )

### set cmake quiet for untyped messages
if    ( NOT SCAI_CMAKE_VERBOSE )
    include ( Functions/setCMakeQuiet )
endif ( NOT SCAI_CMAKE_VERBOSE )

### Set PROJECT NAME / supported LANGUAGES
project ( LAMA_ALL CXX )

## general Settings for all projects
include ( SetBuildFlags ) # VERBOSE, RPATH, LIB64
include ( VersionCheck/Compiler )
include ( Compiler/SetCPPFlags )

## maybe boost
include ( Compiler/CheckC++11 )
if    ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
    include( Package/Boost )
endif ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

include ( Functions/checkValue )
include ( Functions/setAndCheckCache )
include ( Functions/checkWhiteList )
include ( Functions/concatCacheArguments )

include ( VersionDefinition )
include ( SetProjectVariables )
include ( Functions/scaiFunctions )

## Use Boost if BUILD_TEST is enabled
if ( BUILD_TEST )
    include( Package/BoostTest )
endif ( BUILD_TEST )

include ( Dependencies/external )

foreach    ( module ${LAMA_ALL_EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ( module ${LAMA_ALL_EXTERNAL_DEPS} )

include ( Settings/logging )
include ( Settings/tracing )

### finalize build flags
include ( ConcludeFlags )

## summary of what is found
include ( Summaries/scai )

## pass CACHE Variables to all subprojects via white listing

# loads variables WHITELIST, MPI_WHITELIST, CUDA_WHITELIST, MIC_WHITELIST
include ( Dependencies/NeededCacheVariables )

set ( BOOST_ARGS "" )
set ( CUDA_ARGS "" )
set ( GPI2_ARGS "" )
set ( GRAPHPARTITIONING_ARGS "" )
set ( MIC_ARGS "" )
set ( MPI_ARGS "" )
set ( OPENMP_ARGS "" )
set ( SCAI_BLAS_ARGS "" )
set ( JAVA_ARGS "" )

get_cmake_property ( CACHE_VARS CACHE_VARIABLES )
#message ( STATUS "all CACHE_VARS ${CACHE_VARS}")

# always NOTFOUND ?!
#get_cmake_property ( DEFINITION_VARS COMPILE_DEFINITIONS )
#message ( STATUS "all DEFINITION_VARS ${DEFINITION_VARS}")

foreach    ( CACHE_VAR ${CACHE_VARS} )

    #message ( STATUS "test VAR ${CACHE_VAR}" )
    check_whitelist ( CACHE_VAR "${WHITELIST}" "CMAKE_ARGS" )
    
    check_whitelist ( CACHE_VAR "${BOOST_WHITELIST}" "BOOST_ARGS" )
    check_whitelist ( CACHE_VAR "${CUDA_WHITELIST}" "CUDA_ARGS" )
    check_whitelist ( CACHE_VAR "${GPI2_WHITELIST}" "GPI2_ARGS" )
    check_whitelist ( CACHE_VAR "${GRAPHPARTITIONING_WHITELIST}" "GRAPHPARTITIONING_ARGS" )
    check_whitelist ( CACHE_VAR "${MIC_WHITELIST}" "MIC_ARGS" )    
    check_whitelist ( CACHE_VAR "${MPI_WHITELIST}" "MPI_ARGS" )
    check_whitelist ( CACHE_VAR "${GPI_WHITELIST}" "GPI_ARGS" )
    check_whitelist ( CACHE_VAR "${OPENMP_WHITELIST}" "OPENMP_ARGS" )
    check_whitelist ( CACHE_VAR "${SCAI_BLAS_WHITELIST}" "SCAI_BLAS_ARGS" )
    check_whitelist ( CACHE_VAR "${DOXYGEN_WHITELIST}" "DOXYGEN_ARGS" )
    check_whitelist ( CACHE_VAR "${JAVA_WHITELIST}" "JAVA_ARGS" )
    
endforeach ( CACHE_VAR ${CACHE_VARS} )

#message ( STATUS "Passed CMAKE_ARGS ${CMAKE_ARGS}" )

#message ( STATUS "Passed BOOST_ARGS ${BOOST_ARGS}" )
#message ( STATUS "Passed CUDA_ARGS ${CUDA_ARGS}" )
#message ( STATUS "Passed GPI2_ARGS ${GPI_ARGS}" )
#message ( STATUS "Passed GRAPHPARTIONING_ARGS ${GRAPHPARTITIONING_ARGS}" )
#message ( STATUS "Passed MIC_ARGS ${MIC_ARGS}" )
#message ( STATUS "Passed MPI_ARGS ${MPI_ARGS}" )
#message ( STATUS "Passed OPENMP_ARGS ${OPENMP_ARGS}" )
#message ( STATUS "Passed SCAI_BLAS_ARGS ${SCAI_BLAS_ARGS}" )
#message ( STATUS "Passed DOXYGEN_ARGS ${DOXYGEN_ARGS}" )
#message ( STATUS "Passed JAVA_ARGS ${JAVA_ARGS}" )

## define external projects

# loads SCAI_[LIB]_INTERNAL_DEPS variables 
include ( Dependencies/internal )

include ( ExternalProject )

#message ( STATUS "general ${CMAKE_ARGS}" )

## TODO: do in foreach loop
#foreach ( project ${LAMA_ALL_INTERNAL_DEPS} )

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

concatCacheArguments ( scai_hmemo "CACHE_ARGS" )
#message ( STATUS "hmemo CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_hmemo
            DEPENDS ${SCAI_HMEMO_INTERNAL_DEPS}
            PREFIX hmemo
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/hmemo
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR hmemo
            )

concatCacheArguments ( scai_kregistry "CACHE_ARGS" )
#message ( STATUS "kregistry CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_kregistry
            DEPENDS ${SCAI_KREGISTRY_INTERNAL_DEPS}
            PREFIX kregistry
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/kregistry
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR kregistry
            )

concatCacheArguments ( scai_blaskernel "CACHE_ARGS" )
#message ( STATUS "blaskernel CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_blaskernel
            DEPENDS ${SCAI_BLASKERNEL_INTERNAL_DEPS}
            PREFIX blaskernel
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/blaskernel
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR blaskernel
            )
            
concatCacheArguments ( scai_utilskernel "CACHE_ARGS" )
#message ( STATUS "utilskernel CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_utilskernel
            DEPENDS ${SCAI_UTILSKERNEL_INTERNAL_DEPS}
            PREFIX utilskernel
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/utilskernel
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR utilskernel
            )
            
concatCacheArguments ( scai_sparsekernel "CACHE_ARGS" )
#message ( STATUS "sparsekernel CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_sparsekernel
            DEPENDS ${SCAI_SPARSEKERNEL_INTERNAL_DEPS}
            PREFIX sparsekernel
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/sparsekernel
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR sparsekernel
            )

concatCacheArguments ( scai_dmemo "CACHE_ARGS" )
#message ( STATUS "dmemo CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_dmemo
            DEPENDS ${SCAI_DMEMO_INTERNAL_DEPS}
            PREFIX dmemo
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/dmemo
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR dmemo
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

concatCacheArguments ( scai_solver "CACHE_ARGS" )
#message ( STATUS "solver CACHE_ARGS ${CACHE_ARGS}" )

ExternalProject_Add ( scai_solver
            DEPENDS ${SCAI_SOLVER_INTERNAL_DEPS}
            PREFIX solver
            SOURCE_DIR ${CMAKE_SOURCE_DIR}/solver
            CMAKE_ARGS "-DSCAI_COMPLETE_BUILD=1"
            CMAKE_CACHE_ARGS "${CMAKE_ARGS};${CACHE_ARGS}"
            BINARY_DIR solver
            )

if    ( BUILD_DOC )
      add_subdirectory ( doc )
endif ( BUILD_DOC )

## add custom target distclean
include ( CustomCommands/DistClean )

# copy in  scai_test.sh
file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/scai_test.sh
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
       FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )

file ( COPY ${CMAKE_CURRENT_SOURCE_DIR}/scai_test_xml.sh
       DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
       FILE_PERMISSIONS ${PROJECT_FILE_PERMISSIONS} )
