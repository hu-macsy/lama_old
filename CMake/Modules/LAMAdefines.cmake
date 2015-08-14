###
 # @file LAMAdefines.cmake
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
 # @brief Definitions for LAMA
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

## LOGGING Level
#
#  Debug   : use -DLOG_LEVEL_DEBUG
#  Release : use -DLOG_LEVEL_INFO
#  
#  For serious problems: -DLOG_LEVEL_TRACE
#  For benchmarks:       -DLOG_LEVEL_OFF (or -DLOG_LEVEL_FATAL, -DLOG_LEVEL_ERROR)

list ( APPEND LOG_CHOICES "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "OFF" )

if    ( NOT SCAI_LOG_LEVEL )
    if     ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_LOG_LEVEL "INFO" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    else   ( )
        set ( DEFAULT_LOG_LEVEL "TRACE" )
    endif  ( )
endif ( NOT SCAI_LOG_LEVEL )

set ( SCAI_LOG_LEVEL ${DEFAULT_LOG_LEVEL} CACHE STRING
      "Choose level of compile time logging: ${LOG_CHOICES}" )
set ( CACHE SCAI_LOG_LEVEL PROPERTY STRINGS ${LOG_CHOICES} )
checkValue ( ${SCAI_LOG_LEVEL} "${LOG_CHOICES}" )

add_definitions ( -DSCAI_LOG_LEVEL_${SCAI_LOG_LEVEL} )

## ASSERT Level
#
#  Debug   : use -DASSERT_LEVEL_DEBUG
#  Release : use -DASSERT_LEVEL_ERROR
#  
#  For benchmarks:       -DASSERT_LEVEL_OFF

list ( APPEND ASSERT_CHOICES "DEBUG" "ERROR" "OFF" )

if    ( NOT SCAI_ASSERT_LEVEL )
    if     ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_ASSERT_LEVEL "ERROR" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    else   ( )
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    endif  ( )
endif ( NOT SCAI_ASSERT_LEVEL )

set ( SCAI_ASSERT_LEVEL ${DEFAULT_ASSERT_LEVEL} CACHE STRING
      "Choose level of ASSERT: ${ASSERT_CHOICES}" )
set ( CACHE SCAI_ASSERT_LEVEL PROPERTY STRINGS ${ASSERT_CHOICES} )
checkValue ( ${SCAI_ASSERT_LEVEL} "${ASSERT_CHOICES}" )

add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

## LAMA TRACE LEVEL
#
# If TRACE is set to OFF all SCAI_REGION macros in the code are
# completely ignored. If TRACE is set to VT, regions will be traced
# (entry, exit event) for VampirTrace.

list ( APPEND TRACE_CHOICES "ON" "OFF" )

if    ( NOT SCAI_TRACE_LEVEL )
    set ( DEFAULT_TRACE_LEVEL "OFF" )
endif ( NOT SCAI_TRACE_LEVEL )

set ( SCAI_TRACE_LEVEL ${DEFAULT_TRACE_LEVEL} CACHE STRING
     "Choose level of TRACE: ON (enabled) or OFF (disabled, default)" )
set ( CACHE SCAI_TRACE_LEVEL PROPERTY STRINGS ${TRACE_CHOICES} )
checkValue ( ${SCAI_TRACE_LEVEL} "${TRACE_CHOICES}" )

add_definitions ( -DSCAI_TRACE_${SCAI_TRACE_LEVEL} )
