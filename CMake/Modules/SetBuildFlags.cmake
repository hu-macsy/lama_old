###
 # @file Variables.cmake
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief Important CMake variable definitions
 # @author Jan Ecker
 # @date 16.04.2013
 # @since 1.0.0
###

# Check if verbose mode for CMAKE is selected
if    ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    set ( SCAI_FIND_PACKAGE_FLAGS )
else  ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    set ( SCAI_FIND_PACKAGE_FLAGS QUIET )
endif ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )

if    ( SCAI_BUILD_LIB_ONLY )
	set ( BUILD_DOC OFF )
	set ( BUILD_EXAMPLES OFF )
	set ( BUILD_TEST OFF )
endif ( SCAI_BUILD_LIB_ONLY )

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )

## set default switches or check user input

include ( Settings/switchChoices )
include ( Functions/checkValue )

## DOC

# Check if doc should be build
if    ( DEFINED BUILD_DOC )
	#do nothing
else  ( DEFINED BUILD_DOC )
    set ( BUILD_DOC ${BUILD_DOC_DEFAULT} )
endif ( DEFINED BUILD_DOC )
checkValue ( ${BUILD_DOC} "${TRUE_FALSE_CHOICES}" )
message ( "Build doc is set to " ${BUILD_DOC} )

# Choose Doc type
if    ( DEFINED SCAI_DOC_TYPE )
	# do nothing
else  ( DEFINED SCAI_DOC_TYPE )
    set ( SCAI_DOC_TYPE ${SCAI_DOC_TYPE_DEFAULT} )
endif ( DEFINED SCAI_DOC_TYPE )
checkValue ( ${SCAI_DOC_TYPE} "${SCAI_DOC_TYPE_CHOICES}" )
message ( "Doc type is set to " ${SCAI_DOC_TYPE_CHOICES} )

if    ( SCAI_DOC_TYPE STREQUAL json )
    set ( DOC_EXTENTSION "fjson" )
else  ( SCAI_DOC_TYPE STREQUAL json )
    set ( DOC_EXTENTSION "html" )
endif ( SCAI_DOC_TYPE STREQUAL json )

## EXAMPLES

# Check if examples should be build
if    ( DEFINED BUILD_EXAMPLES )
	# do nothing
else  ( DEFINED BUILD_EXAMPLES )
    set ( BUILD_EXAMPLES ${BUILD_EXAMPLES_DEFAULT} )
endif ( DEFINED BUILD_EXAMPLES )
checkValue ( ${BUILD_EXAMPLES} "${TRUE_FALSE_CHOICES}" )
message ( "Build examples is set to " ${BUILD_EXAMPLES} )

## TEST

## Check if tests should be build
if    ( DEFINED BUILD_TEST )
	# do nothing
else  ( DEFINED BUILD_TEST )
    set ( BUILD_TEST ${BUILD_TEST_DEFAULT} )
endif ( DEFINED BUILD_TEST )
checkValue ( ${BUILD_TEST} "${TRUE_FALSE_CHOICES}" )
message ( "Build test is set to " ${BUILD_TEST} )

## CODE COVERAGE

## Check if lama should be build for code coverage
if    ( DEFINED USE_CODE_COVERAGE )
	# do nothing
else  ( DEFINED USE_CODE_COVERAGE )
    set ( USE_CODE_COVERAGE ${USE_CODE_COVERAGE_DEFAULT} )
endif ( DEFINED USE_CODE_COVERAGE )
checkValue ( ${USE_CODE_COVERAGE} "${TRUE_FALSE_CHOICES}" )
message ( "Use code coverage is set to " ${USE_CODE_COVERAGE} )

##  CMAKE_BUILD_TYPE
if    ( DEFINED CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE ) # variable may be defined empty
	# do nothing
else  ( DEFINED CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE )
    set ( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE_DEFAULT} )
endif ( DEFINED CMAKE_BUILD_TYPE AND CMAKE_BUILD_TYPE )
checkValue ( ${CMAKE_BUILD_TYPE} "${CMAKE_BUILD_TYPE_CHOICES}" )
message ( "Build type is set to " ${CMAKE_BUILD_TYPE} )

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

## SCAI_LIBRARY_TYPE ( static or shared )
if    ( DEFINED SCAI_LIBRARY_TYPE )
	# do nothing
else  ( DEFINED SCAI_LIBRARY_TYPE )
    set ( SCAI_LIBRARY_TYPE ${SCAI_LIBRARY_TYPE_DEFAULT} )
endif ( DEFINED SCAI_LIBRARY_TYPE )
checkValue ( ${SCAI_LIBRARY_TYPE} "${SCAI_LIBRARY_TYPE_CHOICES}" )
message ( "Library type is set to " ${SCAI_LIBRARY_TYPE} )
