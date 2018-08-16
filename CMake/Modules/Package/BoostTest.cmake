###
 # @file Package/BoostTest.cmake
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
 # @brief findPackage and configuration of Boost
 # @author Jan Ecker
 # @date 25.04.2013
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

scai_pragma_once ()

### BOOST_INCLUDE_DIR      - Boost include directory
### BOOST_unit_test_framework_FOUND      - if Boost component is found
### Boost_unit_test_framework_LIBRARY    - Boost component library
### BOOST_VERSION          - concluded ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}
### BOOST_TEST_ENABLED     - if Boost_UNIT_TEST_FRAMEWORK_FOUND AND BUILD_TEST

### Boost_USE_STATIC_LIBS  ( default is OFF )
### 
###    - should be set to OFF, default ( will use dynamic Boost libraries if available )
###    - should only be set to ON  if only static BOOST libraries are availabe 
###    - On Linux systems static 64-bit libraries must have been compiled with -fPIC

# set ( Boost_USE_STATIC_LIBS OFF )
# set ( Boost_USE_MULTITHREADED OFF )

set ( BOOST_TEST_MINIMUM_VERSION 1.61 )    ## functionality before is not sufficient

if    ( WIN32 )
    message ( STATUS "Setting special Boost options on Windows" )
    #set ( Boost_USE_STATIC_LIBS ON )
    set ( Boost_USE_MULTITHREADED ON )
endif ( WIN32 )

# Finds packages with custom search options 

set ( Boost_COMPONENTS unit_test_framework )

# FindBoost Debug options comment
if    ( SCAI_CMAKE_VERBOSE )
    set ( Boost_DEBUG TRUE )
    set ( Boost_DETAILED_FAILURE_MSG TRUE )
endif ( SCAI_CMAKE_VERBOSE )

# Find Boost 

find_package ( Boost ${SCAI_FIND_PACKAGE_FLAGS} ${BOOST_TEST_MINIMUM_VERSION} COMPONENTS ${Boost_COMPONENTS} )

set ( BOOST_VERSION "${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" )

if ( Boost_INCLUDE_DIR )
    set ( BOOST_INCLUDE_DIR "${Boost_INCLUDE_DIR}" ) # for getting the module names straight
endif ( Boost_INCLUDE_DIR )

set ( FOUND_BOOST_TEST FALSE )

if ( Boost_UNIT_TEST_FRAMEWORK_FOUND )
    set ( FOUND_BOOST_TEST TRUE )
endif ()

# LAMA irrelevant entries will be removed from cmake GUI completely

set ( Boost_DIR "${Boost_DIR}" CACHE INTERNAL "" )

scai_build_variable ( NAME      USE_BOOST_TEST
                      BOOL 
                      DEFAULT   ${FOUND_BOOST_TEST}
                      DOCSTRING "Boost unit test framework for LAMA unit tests" )

if ( USE_BOOST_TEST )
    if ( NOT FOUND_BOOST_TEST )
        message ( FATAL_ERROR "Forced USE_BOOST_TEST, but not found" )
    endif ()
endif ()

scai_summary_external ( NAME      "Boost Unit Test"
                        ENABLED   ${USE_BOOST_TEST}
                        FOUND     ${FOUND_BOOST_TEST}
                        VERSION   "${BOOST_VERSION}"
                        INCLUDE   "${Boost_INCLUDE_DIR}"
                        LIBRARIES "${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}"
                      )
