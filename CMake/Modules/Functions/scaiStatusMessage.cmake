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

# prints colored text messages
# inspired by soci colormsg function
function ( scai_status_message )
    include ( Settings/bashFormats )

    set ( ERROR "${TextRed}" )
    set ( WARNING "${TextAmber}" )
    set ( INFO "${TextGreen}" )
    set ( HEADLINE "${TextUnderline}" )
    
    set ( coloron FALSE )
    set ( str "" )
    foreach    ( arg ${ARGV} )

        if    ( DEFINED ${arg} )
            if    ( CMAKE_COLOR_MAKEFILE )
                set ( str "${str}${${arg}}" )
                set ( coloron TRUE )
            endif ( CMAKE_COLOR_MAKEFILE )
        else  ( DEFINED ${arg} )
            set ( str "${str}${arg}" )
            if    ( coloron )
                set ( str "${str}${TextColorReset}" )
                set ( coloron FALSE )
            endif ( coloron )
            set ( str "${str} " )
        endif ( DEFINED ${arg} )

    endforeach ( arg ${ARGV} )
    
    message ( STATUS ${str} )
endfunction ( scai_status_message )

function    ( indent_message LEVEL MSG )
    include ( Functions/scaiGenerateBlanks )

    math ( EXPR NUM_BLANKS "2 * (${LEVEL} - 1)" )
    createBlanks ( TYPE_INTENT ${NUM_BLANKS} )

    message ( STATUS "${TYPE_INTENT}${MSG}" )
endfunction ( indent_message LEVEL MSG )

function    ( emptyline )
    message ( STATUS "" )
endfunction ( emptyline )


function    ( format_text )
    # first arg is name of "return variable"
    list ( GET ARGV 0 RESULT_NAME )
    list ( REMOVE_AT ARGV 0 )

    include ( Settings/bashFormats )

    set ( coloron FALSE )
    set ( str "" )
    foreach    ( arg ${ARGV} )

        if    ( DEFINED ${arg} )
            if    ( CMAKE_COLOR_MAKEFILE )
                set ( str "${str}${${arg}}" )
                set ( coloron TRUE )
            endif ( CMAKE_COLOR_MAKEFILE )
        else  ( DEFINED ${arg} )
            set ( str "${str}${arg}" )
            if    ( coloron )
                set ( str "${str}${TextColorReset}" )
                set ( coloron FALSE )
            endif ( coloron )
            set ( str "${str} " )
        endif ( DEFINED ${arg} )

    endforeach ( arg ${ARGV} )

    set ( str "${str}${TextReset}${BGReset}" )
    set ( ${RESULT_NAME} ${str} PARENT_SCOPE )
endfunction ( format_text )

function    ( heading TEXT )
    emptyline()
    format_text ( H1 "TextUnderline" "${TEXT}" )
    indent_message ( "1" "${H1}" )
endfunction ( heading TEXT )

function    ( heading2 TEXT VAR )
    emptyline()

    if    ( VAR STREQUAL "" )
        set ( H2 "" )
    else ( VAR STREQUAL "" )
        if    ( ${VAR} )
            set ( FLAG_TEXT "COMPLETE" )
            set ( FLAG_FORMAT "TextGreen" )
        else  ( ${VAR} )
            set ( FLAG_TEXT "INCOMPLETE" )
            set ( FLAG_FORMAT "TextRed" )
        endif ( ${VAR} )
        format_text ( H2 "${FLAG_FORMAT}" "${FLAG_TEXT}" )
    endif ( VAR STREQUAL "" )
        
    indent_message ( "2" "${TEXT} ${H2}" )
    #emptyline()
endfunction ( heading2 TEXT )

function    ( heading3 TEXT VAR )
    emptyline()
    if    ( ${VAR} )
        set ( FLAG_TEXT "ENABLED" )
        set ( FLAG_FORMAT "TextGreen" )
    else  ( ${VAR} )
        set ( FLAG_TEXT "DISABLED" )
        set ( FLAG_FORMAT "TextAmber" )
    endif ( ${VAR} )
    format_text ( H3 "${FLAG_FORMAT}" "${FLAG_TEXT}" )
        
    indent_message ( "3" "${TEXT} ${H3}" )
endfunction ( heading3 TEXT )

function    ( found_message TEXT VAR OPTIONAL ADDITIONAL_TEXT )
    include ( Functions/scaiGenerateBlanks )

    set ( SCAI_SUMMARY_PACKAGE_NAME_LENGTH 18 )
    scai_generate_blanks ( SCAI_PACKAGE_NAME_BLANKS ${TEXT} ${SCAI_SUMMARY_PACKAGE_NAME_LENGTH} )

    if    ( ${VAR} )
        set ( FLAG_TEXT "FOUND" )
        set ( FLAG_FORMAT "TextGreen" )
    else  ( ${VAR} )
        set ( FLAG_TEXT "NOT FOUND" )
        set ( ADDITIONAL_TEXT "" )
        if     ( ${OPTIONAL} MATCHES "OPTIONAL" )
            set ( FLAG_FORMAT "TextAmber" )
        elseif ( ${OPTIONAL} MATCHES "REQUIRED" )
            set ( FLAG_FORMAT "TextRed" )
        else   ( )
            message ( WARNING "No valid third parameter given to scai_summary_message." )
            set ( FLAG_FORMAT "TextRed" )   
        endif  ( )
    endif ( ${VAR} )
    format_text ( H4 "${FLAG_FORMAT}" "${FLAG_TEXT}" )
        
    indent_message ( "4" "${TEXT}${SCAI_PACKAGE_NAME_BLANKS}${H4} ${ADDITIONAL_TEXT}" )
endfunction ( found_message TEXT )
