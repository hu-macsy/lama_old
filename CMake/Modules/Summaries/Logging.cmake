###
 # @file Summaries/Logging.cmake
 #
 # @license
 # Copyright (c) 2009-2013
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell###
 # @file Summaries/Logging.cmake
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
 # @brief SCAI Logging Configuration Summary
 # @author Lauretta Schubert
 # @date 25.08.2015
 # @since 2.0.0
###

include ( Functions/scaiStatusMessage )
include ( Functions/scaiSummaryMessage )

include ( VersionDefinition )
include ( CompilerVersion )
include ( CheckC++11 )

message ( STATUS "" )
message ( STATUS "Summary of SCAI Logging Configuration:" )
message ( STATUS "=====================================" )
message ( STATUS "" )

scai_status_message ( HEADLINE "Compiler:" )
# C++ Compiler
scai_summary_message ( "FOUND"
                       "CMAKE_CXX_COMPILER"
                       "C++ Compiler"
                       "${CMAKE_CXX_COMPILER_ID} ${${CMAKE_CXX_COMPILER_ID}CXX_COMPILER_VERSION}" )

message ( STATUS "" )

if    ( CXX_SUPPORTS_C11 OR Boost_INCLUDE_DIR )
    set( REQUIRED_FOUND TRUE )
else  ( CXX_SUPPORTS_C11 OR Boost_INCLUDE_DIR )
	set( REQUIRED_FOUND FALSE )
endif ( CXX_SUPPORTS_C11 OR Boost_INCLUDE_DIR )

scai_summary_message ( "STATIC"
                       "REQUIRED_FOUND"
                       "Logging"
                       "Needs compiler supporting C++11 or Boost" )

scai_summary_message ( "FOUND"
					   "CXX_SUPPORTS_C11"
					   "C++11 support"
					   "" )
				
if    ( NOT CXX_SUPPORTS_C11 )
    scai_summary_message ( "FOUND"
                           "Boost_INCLUDE_DIR"
                           "Boost"
                           "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}, add include dir ${Boost_INCLUDE_DIR} to compile your sources" )
endif ( NOT CXX_SUPPORTS_C11 )

# LAMA (core)
message ( STATUS "" )
scai_status_message ( HEADLINE "LIBRARIES:" )

# LAMA CUDA
message ( STATUS "" )
scai_summary_message ( "FOUND"
                       "SCAI_COMMON_FOUND"
                       "SCAI Common"
                       "" )

message ( STATUS "" )

scai_status_message ( HEADLINE "INFO:" )
message ( STATUS "LAMA Version : ${LAMA_VERSION} ${LAMA_VERSION_NAME}" )
message ( STATUS "Build Type   : ${CMAKE_BUILD_TYPE}" )
message ( STATUS "Library Type : ${SCAI_LIBRARY_TYPE}" )
message ( STATUS "ASSERT Level : ${SCAI_ASSERT_LEVEL}" )
message ( STATUS "" )
