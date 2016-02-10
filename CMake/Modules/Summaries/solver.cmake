###
 # @file Summaries/solver.cmake
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
 # @brief LAMA Summary for build configuration
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

include ( Functions/scaiStatusMessage )
include ( Functions/scaiSummaryMessage )

message ( STATUS "" )
message ( STATUS "Summary of solver Configuration:" )
message ( STATUS "==============================" )
message ( STATUS "" )

scai_status_message ( HEADLINE "Compiler:" )
# C++ Compiler
scai_summary_message ( "FOUND"
                       "CMAKE_CXX_COMPILER"
                       "C++ Compiler"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )
                       
message ( STATUS "" )

if    ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR )
    set( REQUIRED_FOUND TRUE )
else  ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR )
	set( REQUIRED_FOUND FALSE )
endif ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "solver"
                       "Needs compiler supporting C++11 or Boost and pThreads" )

scai_summary_message ( "FOUND"
					             "CXX_SUPPORTS_C11"
					             "C++11 support"
					             "" )

if    ( NOT CXX_SUPPORTS_C11 )
    scai_summary_message ( "FOUND"
                           "BOOST_INCLUDE_DIR"
                           "Boost"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}, add include dir ${BOOST_INCLUDE_DIR} to compile your sources" )
endif ( NOT CXX_SUPPORTS_C11 )

message ( STATUS "" )

# LAMA (core)
message ( STATUS "" )
scai_status_message ( HEADLINE "LIBRARIES:" )
   
set ( REQUIRED_FOUND FALSE )
if    ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_HMEMO_FOUND AND 
        SCAI_KREGISTRY_FOUND AND SCAI_BLASKERNEL_FOUND AND SCAI_LAMA_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( )

message ( STATUS "" )
scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "Internal Libraries (core)"
                       "" )

    scai_summary_message ( "FOUND"
                           "SCAI_COMMON_FOUND"
                           "SCAI common"
                           "" )
                           
    scai_summary_message ( "FOUND"
                           "SCAI_LOGGING_FOUND"
                           "SCAI logging"
                           "" )
                           
    scai_summary_message ( "FOUND"
                           "SCAI_TRACING_FOUND"
                           "SCAI tracing"
                           "" )

    scai_summary_message ( "FOUND"
                           "SCAI_TASKING_FOUND"
                           "SCAI tasking"
                           "" )
                               
    scai_summary_message ( "FOUND"
                           "SCAI_HMEMO_FOUND"
                           "SCAI hmemo"
                           "" )
                           
    scai_summary_message ( "FOUND"
                           "SCAI_KREGISTRY_FOUND"
                           "SCAI kregistry"
                           "" )

    scai_summary_message ( "FOUND"
                           "SCAI_BLASKERNEL_FOUND"
                           "SCAI blaskernel"
                           "" )

    scai_summary_message ( "FOUND"
                           "SCAI_LAMA_FOUND"
                           "SCAI lama"
                           "" )

# LAMA TEST
scai_status_message ( HEADLINE "TESTING:" )

scai_summary_message ( "USE"
                       "BUILD_TEST"
                       "LAMA TEST"
                       "" )

    # Boost Test-Framework
    scai_summary_message ( "FOUND"
                           "Boost_UNIT_TEST_FRAMEWORK_FOUND"
                           "Boost Unit Test"
                           "" )
                           
    # Boost Regex
    scai_summary_message ( "FOUND"
                           "Boost_REGEX_FOUND"
                           "Boost Regex"
                           "" )

message ( STATUS "" )

scai_status_message ( HEADLINE "INFO:" )

message ( STATUS "Solver Version : ${SCAI_SOLVER_VERSION} ${SCAI_VERSION_NAME}" )
message ( STATUS "Build Type   : ${CMAKE_BUILD_TYPE}" )
message ( STATUS "Library Type : ${SCAI_LIBRARY_TYPE}" )
message ( STATUS "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
message ( STATUS "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" )
message ( STATUS "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" )
if    ( USE_CODE_COVERAGE )
	message ( STATUS "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )
message ( STATUS "" )
