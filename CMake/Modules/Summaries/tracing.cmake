###
 # @file Summaries/tracing.cmake
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
 # @brief SCAI tracing Configuration Summary
 # @author Lauretta Schubert
 # @date 25.08.2015
 # @since 2.0.0
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "=====================================" )
message ( STATUS "Summary of SCAI tracing Configuration:" )
message ( STATUS "=====================================" )

include ( Summaries/Modules/Compiler )

# tracing (core)
set ( REQUIRED_FOUND FALSE )
if    ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND )

heading2 ( "Required core" "REQUIRED_FOUND" )
heading3 ( "Internal Libraries" "REQUIRED_FOUND" )
    found_message ( "SCAI common"       "SCAI_COMMON_FOUND"       "REQUIRED" "Version ${SCAI_COMMON_VERSION}"       )
    found_message ( "SCAI logging"      "SCAI_LOGGING_FOUND"      "REQUIRED" "Version ${SCAI_LOGGING_VERSION}"      )

heading2 ( "Optional components" "" )
heading3 ( "Gui:" "Java_FOUND" )
    found_message ( "Java" "Java_FOUND" "OPTIONAL" "" )

include ( Summaries/Modules/Configuration )