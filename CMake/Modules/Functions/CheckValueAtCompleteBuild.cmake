###
 # @file Functions.cmake
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
 # @brief CMake function for option checking if SCAI_COMPLETE_BUILD 
 # @author Lauretta Schubert
 # @date 28.09.2015
 # @since 2.0.0
###
include ( Functions/checkValue )

function ( checkValueAtCompleteBuild LIBRARY )

string ( TOUPPER ${LIBRARY} UPPER_LIBRARY )
set ( VAR_NAME "SCAI_${UPPER_LIBRARY}_EXTERNAL_DEPS" )

## from SetBuildFlags

# CMAKE_BUILD_TYPE
set ( CMAKE_BUILD_TYPE_OPTIONS None Debug Release RelWithDebInfo MinSizeRel )
checkValue ( ${CMAKE_BUILD_TYPE} "${CMAKE_BUILD_TYPE_OPTIONS}" )

# SCAI_LIBRARY_TYPE
set ( SCAI_LIBRARY_TYPE_OPTIONS STATIC SHARED )
checkValue ( ${SCAI_LIBRARY_TYPE} "${SCAI_LIBRARY_TYPE_OPTIONS}" )

## from SCAIAssert
list ( APPEND ASSERT_CHOICES "DEBUG" "ERROR" "OFF" )
checkValue ( ${SCAI_ASSERT_LEVEL} "${ASSERT_CHOICES}" )

## from SCAI_BLAS
list ( FIND ${VAR_NAME} SCAI_BLAS test_var )
message ( STATUS "test_var SCAI_BLAS ${test_var}" )
if    ( ${test_var} GREATER -1 )
	LIST ( APPEND LIBRARY_CHOICES "auto" "MKL" "BLAS" "INTERNALBLAS" )
	checkValue( ${SCAI_BLAS_LIBRARY} "${LIBRARY_CHOICES}" )
endif ( ${test_var} GREATER -1 )

## from SetNVCCFlags
list ( FIND ${VAR_NAME} CUDA test_var )
message ( STATUS "test_var CUDA ${test_var}" )
if    ( ${test_var} GREATER -1 )
	list ( APPEND CC_CHOICES "not-found" "13" "20" "21" "30" "32" "35" "50" )
	checkValue( ${CUDA_COMPUTE_CAPABILITY} "${CC_CHOICES}" )
endif ( ${test_var} GREATER -1 )

endfunction ( checkValueAtCompleteBuild LIBRARY )