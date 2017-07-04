###
 # @file SetBuildFlags.cmake
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
 # @brief Important CMake variable definitions
 # @author Thomas Brandes, Lauretta Schubert, Jan Ecker
 # @date 16.04.2013
###

# Check if verbose mode for CMAKE is selected
if ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    set ( SCAI_FIND_PACKAGE_FLAGS )
else ()
    set ( SCAI_FIND_PACKAGE_FLAGS QUIET )
endif ()

if ( SCAI_BUILD_LIB_ONLY )
    set ( BUILD_DOC OFF )
    set ( BUILD_EXAMPLES OFF )
    set ( BUILD_TEST OFF )
endif ( SCAI_BUILD_LIB_ONLY )

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )

## set default switches or check user input

include ( Settings/switchChoices )
include ( Functions/checkValue )
include ( Functions/listToString )
include ( Functions/parseBoolean )

## DOC

include( Package/doc )

# Check if doc should be build
if    ( DEFINED BUILD_DOC )
    parseBoolean( BUILD_DOC )

    if    ( BUILD_DOC AND NOT DOC_FOUND )
        message( FATAL_ERROR "Build of documentation enabled, but configuration is incomplete!")
    endif ( BUILD_DOC AND NOT DOC_FOUND )

else  ( DEFINED BUILD_DOC )
    
    if    ( DOC_FOUND )
        set ( BUILD_DOC ON )
    else  ( DOC_FOUND )
        set ( BUILD_DOC OFF )
    endif ( DOC_FOUND )

endif ( DEFINED BUILD_DOC )
checkValue ( ${BUILD_DOC} "${TRUE_FALSE_CHOICES}" )
set ( BUILD_DOC ${BUILD_DOC} CACHE BOOL "Enable / Disable building of doc" )

set( DOC_ENABLED OFF )
if     ( DOC_FOUND AND BUILD_DOC )
    set( DOC_ENABLED ON )
endif ( DOC_FOUND AND BUILD_DOC )

# Choose Doc type
if    ( DEFINED SCAI_DOC_TYPE )
    # do nothing
else  ( DEFINED SCAI_DOC_TYPE )
    set ( SCAI_DOC_TYPE ${SCAI_DOC_TYPE_DEFAULT} )
endif ( DEFINED SCAI_DOC_TYPE )
checkValue ( ${SCAI_DOC_TYPE} "${SCAI_DOC_TYPE_CHOICES}" )
set ( SCAI_DOC_TYPE ${SCAI_DOC_TYPE} CACHE STRING "Choose the type of documentation, options are: ${SCAI_DOC_TYPE_CHOICES}." )

## EXAMPLES

# Check if examples should be build
if    ( DEFINED BUILD_EXAMPLES )
    parseBoolean( BUILD_EXAMPLES )
else  ( DEFINED BUILD_EXAMPLES )
    set ( BUILD_EXAMPLES ${BUILD_EXAMPLES_DEFAULT} )
endif ( DEFINED BUILD_EXAMPLES )

checkValue ( ${BUILD_EXAMPLES} "${TRUE_FALSE_CHOICES}" )

set ( BUILD_EXAMPLES ${BUILD_EXAMPLES} CACHE BOOL "Enable / Disable building of examples" )

## TEST

## Check if tests should be build
if    ( DEFINED BUILD_TEST )
    parseBoolean( BUILD_TEST )
else  ( DEFINED BUILD_TEST )
    set ( BUILD_TEST ${BUILD_TEST_DEFAULT} )
endif ( DEFINED BUILD_TEST )
checkValue ( ${BUILD_TEST} "${TRUE_FALSE_CHOICES}" )
set ( BUILD_TEST ${BUILD_TEST} CACHE BOOL "Enable / Disable building of tests" )

## CODE COVERAGE

## Check if lama should be build for code coverage
if ( DEFINED USE_CODE_COVERAGE )
    # do nothing
    parseBoolean( USE_CODE_COVERAGE )
else  ( DEFINED USE_CODE_COVERAGE )
    set ( USE_CODE_COVERAGE ${USE_CODE_COVERAGE_DEFAULT} )
endif ( DEFINED USE_CODE_COVERAGE )
checkValue ( ${USE_CODE_COVERAGE} "${TRUE_FALSE_CHOICES}" )
set ( USE_CODE_COVERAGE ${USE_CODE_COVERAGE} CACHE BOOL "Enable / Disable use of Code Coverage" )

##  CMAKE_BUILD_TYPE
if    ( DEFINED CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE ) # variable may be defined empty
    # do nothing
else  ( DEFINED CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE_DEFAULT} )
endif ( DEFINED CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE )
checkValue ( ${CMAKE_BUILD_TYPE} "${CMAKE_BUILD_TYPE_CHOICES}" )
set ( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: ${CMAKE_BUILD_TYPE_CHOICES}." )

## SCAI_ASSERT_LEVEL
if    ( NOT SCAI_ASSERT_LEVEL )
    if     ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( SCAI_ASSERT_LEVEL "ERROR" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( SCAI_ASSERT_LEVEL "DEBUG" )
    else   ( )
        set ( SCAI_ASSERT_LEVEL "DEBUG" )
    endif  ( )
endif ( NOT SCAI_ASSERT_LEVEL )
checkValue ( ${SCAI_ASSERT_LEVEL} "${SCAI_ASSERT_CHOICES}" )
set ( SCAI_ASSERT_LEVEL ${SCAI_ASSERT_LEVEL} CACHE STRING "Choose level of ASSERT: ${SCAI_ASSERT_CHOICES}" )

## SCAI_LIBRARY_TYPE ( static or shared )
if    ( DEFINED SCAI_LIBRARY_TYPE )
    # do nothing
else  ( DEFINED SCAI_LIBRARY_TYPE )
    set ( SCAI_LIBRARY_TYPE ${SCAI_LIBRARY_TYPE_DEFAULT} )
endif ( DEFINED SCAI_LIBRARY_TYPE )
checkValue ( ${SCAI_LIBRARY_TYPE} "${SCAI_LIBRARY_TYPE_CHOICES}" )
set ( SCAI_LIBRARY_TYPE ${SCAI_LIBRARY_TYPE} CACHE STRING "Choose the type of linking: ${SCAI_LIBRARY_TYPE_CHOICES}" )
