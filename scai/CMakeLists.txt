###
 # @file CMakeLists.txt.tb
 #
 # @license
 # Copyright (c) 2009-2017
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
 # @brief CMake configuration file for involved SCAI projects.
 # @author Thomas Brandes
 # @date 03.07.2017
###

cmake_minimum_required ( VERSION 2.8.8 )

project ( LAMA_ALL CXX )

##  Define the list of all SCAI module projects
##  NOTE: list can be extended

set ( SCAI_ALL_MODULES common logging tracing tasking kregistry hmemo 
                       blaskernel utilskernel sparsekernel dmemo lama solver  )
 
### Define modules that should be used ( inclusive all dependent ones )
### Default: use all modules; but it can be changed via cmake gui
### dependencies are not resolved here but will be verified during configuration

if ( DEFINED SCAI_MODULES )
    set ( SCAI_USED_MODULES ${SCAI_MODULES} CACHE STRING "Used SCAI module projects" )
else ()
    set ( SCAI_USED_MODULES ${SCAI_ALL_MODULES} CACHE STRING "Used SCAI module projects" )
endif ()

### Define some general variables used for identifaction 

set ( SCAI_VERSION_NAME   "Barrancas Blancas" )
set ( SCAI_COPYRIGHT_YEAR 2017 )
set ( SCAI_VERSION        2.2.0 )
set ( SCAI_LIBRARY_PREFIX "scai" )             # more unique names for all SCAI module libraries

### Set Module Path

set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../CMake/Modules" )

### set cmake quiet for untyped messages
if    ( NOT SCAI_CMAKE_VERBOSE )
    include ( scai_function/setCMakeQuiet )
endif ( NOT SCAI_CMAKE_VERBOSE )

### General stuff for all projects

include ( VersionCheck/Compiler )
include ( Compiler/SetCPPFlags )

include ( scai_macro/scai_summary )

scai_summary_external ( NAME       CXX
                        ENABLED    True
                        FOUND      True
                        VERSION    "${CMAKE_CXX_COMPILER_ID} ${CXX_COMPILER_VERSION}"
                        EXECUTABLE "${CMAKE_CXX_COMPILER}"
                        CXX_FLAGS  "${CMAKE_CXX_COMPILER_FLAGS}"
                      )

## C++11 or Boost

include ( Compiler/CheckC++11 )

if ( CXX_SUPPORTS_C11 )

    scai_summary_external ( NAME      "C++11 support"
                            ENABLED   True
                            FOUND     True
                            VERSION   "via C++ compiler"
                          )
else ()

    include( Package/Boost )

    scai_summary_external ( NAME      "C++11 Support"
                            ENABLED   True
                            FOUND     True
                            VERSION   "via Boost ${BOOST_VERSION}"
                            INCLUDE   "${BOOST_INCLUDE_DIR}"
                          )

    include_directories( ${BOOST_INCLUDE_DIR} )

endif ()

## Define all variables (in CACHE) that are used for individual choices 

include ( scai_macro/scai_build_variable )

scai_build_variable ( NAME      CMAKE_BUILD_TYPE
                      CHOICES   "Debug" "Release" "RelWithDebInfo" "MinSizeRel" ""
                      DEFAULT   "Debug"
                      DOCSTRING "Set build type" )

scai_build_variable ( NAME      USE_CODE_COVERAGE
                      BOOL
                      DEFAULT   False
                      DOCSTRING "use of code coverage" )

scai_build_variable ( NAME      BUILD_TEST
                      BOOL
                      DEFAULT   True
                      DOCSTRING "build of test" )

scai_build_variable ( NAME      BUILD_EXAMPLES
                      BOOL
                      DEFAULT   True
                      DOCSTRING "build of examples" )

scai_build_variable ( NAME      BUILD_DOC
                      BOOL
                      DEFAULT   True
                      DOCSTRING "build of user doc" )

scai_build_variable ( NAME      SCAI_LIBRARY_TYPE
                      CHOICES   STATIC SHARED
                      DEFAULT   SHARED
                      DOCSTRING "Choose the type of linking" )

scai_build_variable ( NAME      SCAI_DOC_TYPE 
                      CHOICES   "html" "json" "xml" "latex"
                      DEFAULT   html
                      DOCSTRING "Choose type of user documentation" )

if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
    set ( SCAI_ASSERT_DEFAULT ERROR )
else ()
    set ( SCAI_ASSERT_DEFAULT DEBUG )
endif ()

scai_build_variable ( NAME      SCAI_ASSERT_LEVEL
                      CHOICES   DEBUG ERROR OFF
                      DEFAULT   ${SCAI_ASSERT_DEFAULT}
                      DOCSTRING "Choose minimal level of assert @ compile time" )
   
if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
    set ( SCAI_LOGGING_DEFAULT "INFO" )
else ()
    set ( SCAI_LOGGING_DEFAULT "DEBUG" )
endif ()

scai_build_variable ( NAME      SCAI_LOGGING_LEVEL
                      CHOICES   "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "OFF"
                      DEFAULT   ${SCAI_LOGGING_DEFAULT}
                      DOCSTRING "Choose minimal level of logging @ compile time" )

scai_build_variable ( NAME      SCAI_TRACING
                      BOOL 
                      DEFAULT   On
                      DOCSTRING "tracing of regions for performance analysis" )

## BUILD_TEST is enabled -> search BoostTest, load macro files

enable_testing ()

## make check : builds all tests and runs them

add_custom_target ( check COMMAND ${CMAKE_CTEST_COMMAND} )

add_custom_target ( examples )

if ( BUILD_TEST )
    include( Package/BoostTest )
    include( scai_macro/scai_test_scripts )  # macro used in test directories
    include( scai_macro/scai_unit_test    )  # macro used in test directories
endif ()

if ( BUILD_EXAMPLES )
    include( scai_macro/scai_run_script        )  # macro used in example directories
    include( scai_macro/scai_example_makefile  )  # macro used in example directories
endif ()

### Sphinx is used for user documentation and a module project 
### does not have to specify it as external dependency

include ( Package/Sphinx )

### CUDA and OpenMP are required to set compiler flags

include ( Package/OpenMP )
include ( Package/CUDA )

include ( ConcludeFlags )

### Traverse all module projects for configuration

foreach ( module ${SCAI_USED_MODULES} )

    # Note: names of subdirectories are same as module names

    add_subdirectory( ${module} )

endforeach ()

include ( Package/Doxygen )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

set ( MODULE_NAME "lama_all" )
set ( INTERNAL_DEPS ${SCAI_USED_MODULES} )

include ( CustomCommands/SphinxDoc )

## add custom target doxygendoc

include ( CustomCommands/DoxygenDoc )

##  add custom target doc to build both

add_custom_target ( doc DEPENDS doxygendoc doc_lama_all )

## add custom target distClean

include ( CustomCommands/DistClean )

## summary of what is found

include ( Summary )

