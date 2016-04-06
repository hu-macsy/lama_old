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
 # @brief CMake functions and macros
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

include ( Functions/scaiGenerateBlanks )

## Need to be macros not functions, because modifications of the parent scope

# generates messages for scai summary page
macro    ( scai_summary_message MESSAGE_TYPE EXPRESSION PACKAGE_NAME ADDITIONAL_INFO )
    set ( SCAI_SUMMARY_PACKAGE_NAME_LENGTH 18 )
    if    ( ${MESSAGE_TYPE} STREQUAL "FOUND" )
        set ( TYPE_TRUE "FOUND" )
        set ( TYPE_FALSE "NOT FOUND" )
        set ( TYPE_INTENT "    " )
        scai_generate_blanks ( SCAI_PACKAGE_NAME_BLANKS ${PACKAGE_NAME} ${SCAI_SUMMARY_PACKAGE_NAME_LENGTH} )
    endif ( ${MESSAGE_TYPE} STREQUAL "FOUND" )
    
    if ( ${MESSAGE_TYPE} STREQUAL "STATIC" )
        set ( TYPE_TRUE "REQUIRED" )
        set ( TYPE_FALSE "REQUIRED" )
        set ( TYPE_INTENT " " )
        set ( SCAI_PACKAGE_NAME_BLANKS "" )
    endif ( ${MESSAGE_TYPE} STREQUAL "STATIC" )
    
    if ( ${MESSAGE_TYPE} STREQUAL "USE" )
        set ( TYPE_TRUE "ENABLED" )
        set ( TYPE_FALSE "DISABLED" )
        set ( TYPE_INTENT "  " )
        set ( SCAI_PACKAGE_NAME_BLANKS "" )
    endif ( ${MESSAGE_TYPE} STREQUAL "USE" )

    if ( ${MESSAGE_TYPE} STREQUAL "HEADLINE" )
        set ( TYPE_TRUE "OK" )
        set ( TYPE_FALSE "FAILED" )
        set ( TYPE_INTENT "" )
        set ( SCAI_PACKAGE_NAME_BLANKS "" )
    endif ( ${MESSAGE_TYPE} STREQUAL "HEADLINE" )

#    if ( DEFINED ${EXPRESSION} )

        if    ( ${EXPRESSION} )
            scai_status_message ( ${TYPE_INTENT} ${PACKAGE_NAME} ${SCAI_PACKAGE_NAME_BLANKS} INFO ${TYPE_TRUE} ${ADDITIONAL_INFO} )
        else  ( ${EXPRESSION} )
            scai_status_message ( ${TYPE_INTENT} ${PACKAGE_NAME} ${SCAI_PACKAGE_NAME_BLANKS} ERROR ${TYPE_FALSE} )  
        endif ( ${EXPRESSION} )

#    else  ( DEFINED ${EXPRESSION} )

#        scai_status_message ( ${TYPE_INTENT} ${PACKAGE_NAME} ${SCAI_PACKAGE_NAME_BLANKS} INFO ${TYPE_TRUE} )

#    endif ( DEFINED ${EXPRESSION} )
endmacro ( scai_summary_message )