###
 # @file scai/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief CMake main configuration file for SCAI LAMA with all modules.
 # @author Thomas Brandes
 # @date 03.07.2017
###

cmake_minimum_required ( VERSION 2.8.8 )

project ( LAMA_ALL CXX )

##  Define the list of all SCAI LAMA modules (subprojects)
##  (should match the list of all subdirectories)
##
##  NOTE: list must be extended if new modules are added

set ( SCAI_ALL_MODULES common logging tracing tasking kregistry hmemo 
                       blaskernel utilskernel sparsekernel dmemo lama 
                       partitioning solver tools )
 
### Define modules that should be used; default: use all modules
### The used modules can be configured explicitly via command line or cmake GUI

if ( DEFINED SCAI_MODULES )
    set ( SCAI_USED_MODULES ${SCAI_MODULES} CACHE STRING "Used SCAI module projects" )
else ()
    set ( SCAI_USED_MODULES ${SCAI_ALL_MODULES} CACHE STRING "Used SCAI module projects" )
endif ()

### Define some general variables used for identifaction 

set ( SCAI_VERSION_NAME   "Chimborazo" )
set ( SCAI_COPYRIGHT_YEAR 2018 )
set ( SCAI_VERSION        3.0.0 )
set ( SCAI_LIBRARY_PREFIX "scai" )             # more unique names for all SCAI module libraries

### Extend CMake module path for SCAI/LAMA specific extensions

set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../CMake/Modules" )

### set cmake quiet for untyped messages
if    ( NOT SCAI_CMAKE_VERBOSE )
    include ( scai_function/setCMakeQuiet )
endif ( NOT SCAI_CMAKE_VERBOSE )

### Setting cmake policies to ignore warnings

# set this policy (from version 3.1.12) 
# cmake_policy ( SET CMP0074 NEW )

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

    message ( FATAL_ERROR "Compiler does not support C++11, which is a necessary requirement." )

endif ()

## Define all variables (in CACHE) that are used for individual choices 

include ( scai_macro/scai_build_variable )

scai_build_variable ( NAME      CMAKE_BUILD_TYPE
                      CHOICES   "Debug" "Release" "RelWithDebInfo" "MinSizeRel" ""
                      DEFAULT   "Release"
                      DOCSTRING "Set build type" )

scai_build_variable ( NAME      USE_CODE_COVERAGE
                      BOOL
                      DEFAULT   False
                      DOCSTRING "use of code coverage" )

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
    set ( SCAI_LOG_DEFAULT "INFO" )
else ()
    set ( SCAI_LOG_DEFAULT "DEBUG" )
endif ()

scai_build_variable ( NAME      SCAI_LOG_LEVEL
                      CHOICES   "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "OFF"
                      DEFAULT   ${SCAI_LOG_DEFAULT}
                      DOCSTRING "Choose minimal level of logging @ compile time" )

scai_build_variable ( NAME      SCAI_TRACE
                      BOOL 
                      DEFAULT   On
                      DOCSTRING "tracing of regions for performance analysis" )

enable_testing ()

## make check : builds all tests and runs them

add_custom_target ( tests )
add_custom_target ( check COMMAND ${CMAKE_CTEST_COMMAND} DEPENDS tests )
add_custom_target ( examples )

## include test macros, packages here so it must not be redone in each test directory

include( Package/BoostTest )
include( scai_macro/scai_test_scripts )  # macro used in test directories
include( scai_macro/scai_test         )  # macro used in test directories

#
file (COPY ${CMAKE_SOURCE_DIR}/run_tests.py  DESTINATION ${CMAKE_BINARY_DIR} )

### Sphinx is used for user documentation, individually for all modules + main
### Doxygen is used for system documentation generated by all sources

include ( Package/Sphinx )
include ( Package/Doxygen )

### CUDA and OpenMP are required to set compiler flags

include ( Package/OpenMP )
include ( Package/CUDA )

include ( ConcludeFlags )

add_subdirectory("testsupport" EXCLUDE_FROM_ALL)

### Traverse all module projects for configuration

foreach ( module ${SCAI_USED_MODULES} )

    # Note: names of subdirectories are same as module names

    add_subdirectory( ${module} )

endforeach ()

## Configure cmake config file SCAIConfig.cmake and install it
## Installation directory share/cmake/SCAI is one installation prefix
## that will be used by CMake when find_package( SCAI .... ) is called

configure_file ( SCAIConfig.cmake.in SCAIConfig.cmake @ONLY )
install ( FILES ${CMAKE_BINARY_DIR}/SCAIConfig.cmake DESTINATION share/cmake/SCAI )

##############################################################################
#                                                                            #
#  Documentation                                                             #
#                                                                            #
##############################################################################

file (COPY ${CMAKE_SOURCE_DIR}/doc/LAMA.html DESTINATION ${CMAKE_BINARY_DIR}/doc )
file (COPY ${CMAKE_SOURCE_DIR}/doc/lama.png  DESTINATION ${CMAKE_BINARY_DIR}/doc )
file (COPY ${CMAKE_SOURCE_DIR}/doc/lgplv3-147x51.png DESTINATION ${CMAKE_BINARY_DIR}/doc )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
## depends on user doc of all defined modules ( = used modules + dependent ones)

set ( MODULE_NAME "libama" )
set ( INTERNAL_DEPS ${SCAI_DEFINED_MODULES} )

## MODULE_NAME, INTERNAL_DEPS are needed to build LAMA user documentation

include ( CustomCommands/SphinxDoc )

set ( DOC_BINARY_DIR ${CMAKE_BINARY_DIR}/doc/user )
set ( DOC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/share/doc/user )

string ( REGEX REPLACE "/" "\\\\/" DOC_BINARY_DIR1  ${DOC_BINARY_DIR} )
string ( REGEX REPLACE "/" "\\\\/" DOC_INSTALL_DIR1 ${DOC_INSTALL_DIR} )

configure_file( doc/rename.sh.in doc/user/rename.sh @ONLY )

install( DIRECTORY ${DOC_BINARY_DIR} DESTINATION share/doc )

## add custom target doxygendoc

include ( CustomCommands/DoxygenDoc )

##  add custom target doc to build both

add_custom_target ( doc DEPENDS doxygendoc doc_libama )

## add custom target distClean

include ( CustomCommands/DistClean )

## summary of what is found

include ( Summary )

