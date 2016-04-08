###
 # @file Summaries/common.cmake
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
 # @brief SCAI common Configuration Summary
 # @author Lauretta Schubert
 # @date 25.08.2015
 # @since 2.0.0
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "=====================================" )
message ( STATUS "Summary of SCAI common Configuration:" )
message ( STATUS "=====================================" )

include ( Summaries/Modules/Compiler )

# common (core)
heading ( "External Libraries:" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR )
    set ( REQUIRED_FOUND TRUE )
endif ( SCAI_THREAD_LIBRARIES AND SCAI_BOOST_INCLUDE_DIR  )

heading2 ( "Required core" "REQUIRED_FOUND" )

    # pthreads
    found_message ( "pThreads" "SCAI_THREAD_LIBRARIES" "REQUIRED" "Version ${SCAI_THREAD_VERSION}" )
    # boost
    found_message ( "Boost" "SCAI_BOOST_INCLUDE_DIR" "REQUIRED" "Version ${BOOST_VERSION} at ${SCAI_BOOST_INCLUDE_DIR}" )

heading2 ( "Optional components" "" )
include ( Summaries/Modules/Build )

#include ( Summaries/Modules/Configuration )
heading ( "Configuration Details:" )

indent_message ( "1" "SCAI common Version : ${SCAI_COMMON_VERSION}" )
indent_message ( "1" "Build Type   : ${CMAKE_BUILD_TYPE}" )
indent_message ( "1" "Library Type : ${SCAI_LIBRARY_TYPE}" )
indent_message ( "1" "ASSERT Level : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
#indent_message ( "1" "LOG Level    : ${SCAI_LOGGING_LEVEL} ( -D${SCAI_LOGGING_FLAG} )" ) #opt
#indent_message ( "1" "TRACING      : ${SCAI_TRACING} ( -D${SCAI_TRACING_FLAG} )" ) #opt
if    ( USE_CODE_COVERAGE )
    indent_message ( "1" "CODE COVERAGE: ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )

emptyline()
