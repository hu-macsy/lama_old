###
 # @file scaiMessages.cmake
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
 # @brief CMake functions for formatted messages used in summaries
 # @author Jan Ecker
 # @date 25.04.2013
###

##  Help function 

function    ( createBlanks OUTPUT LENGTH )
    set ( MESSAGE_BLANKS "")
    foreach    ( SCAI_I RANGE ${LENGTH} )
        set ( MESSAGE_BLANKS "${MESSAGE_BLANKS} " )
    endforeach ( SCAI_I RANGE ${LENGTH} )

    set ( ${OUTPUT} ${MESSAGE_BLANKS} PARENT_SCOPE )
endfunction ( createBlanks OUTPUT LENGTH )

function    ( scai_generate_blanks OUTPUT STRING LENGTH )
    string ( LENGTH "${STRING}" SCAI_STRING_LENGTH )
    # -1 for correct looping from 0 to LENGTH
    math ( EXPR SCAI_MESSAGE_BLANK_LENGTH ${LENGTH}-${SCAI_STRING_LENGTH} )
    
    createBlanks ( MESSAGE_BLANKS ${SCAI_MESSAGE_BLANK_LENGTH} )
    
    set ( ${OUTPUT} ${MESSAGE_BLANKS} PARENT_SCOPE )
endfunction ( scai_generate_blanks )

## Set all escape sequences for bash formatting

if    ( NOT WIN32 )

    # full bash color and formating escape sequences:
    # http://misc.flogisoft.com/bash/tip_colors_and_formatting

    string ( ASCII 27 Esc )

    set ( TextReset          "${Esc}[0m" )

    set ( TextBold           "${Esc}[1m" )
    set ( TextDim            "${Esc}[2m" )
    set ( TextItalic         "${Esc}[3m" )
    set ( TextUnderline      "${Esc}[4m" )
    set ( TextBlink          "${Esc}[5m" )
    set ( TextReverse        "${Esc}[7m" )
    set ( TextHidden         "${Esc}[8m" )

    set ( TextBoldReset      "${Esc}[21m" )
    set ( TextDimReset       "${Esc}[22m" )
    set ( TextItalicReset    "${Esc}[23m")
    set ( TextUnderlineReset "${Esc}[24m" )
    set ( TextBlinkReset     "${Esc}[25m" )
    set ( TextReverseReset   "${Esc}[27m" )
    set ( TextHiddenReset    "${Esc}[28m" )

    set ( TextColorReset     "${Esc}[m" )
    set ( TextBlack          "${Esc}[30m" )
    set ( TextRed            "${Esc}[31m" )
    set ( TextGreen          "${Esc}[32m" )
    set ( TextYellow         "${Esc}[33m" )
    set ( TextBlue           "${Esc}[34m" )
    set ( TextMagenta        "${Esc}[35m" )
    set ( TextCyan           "${Esc}[36m" )  
    set ( TextLightGray      "${Esc}[37m" )
    set ( TextDarkGray       "${Esc}[90m" )
    set ( TextLightRed       "${Esc}[91m" )
    set ( TextLightGreen     "${Esc}[92m" )
    set ( TextLightYellow    "${Esc}[93m" )
    set ( TextLightBlue      "${Esc}[94m" )
    set ( TextLightMagenta   "${Esc}[95m" )
    set ( TextLightCyan      "${Esc}[96m" )
    set ( TextWhite          "${Esc}[97m" )

    set ( BGColorReset       "${Esc}[49m" )

    set ( BGBlack            "${Esc}[40m" )
    set ( BGRed              "${Esc}[41m" )
    set ( BGGreen            "${Esc}[42m" )
    set ( BGYellow           "${Esc}[43m" )
    set ( BGBlue             "${Esc}[44m" )
    set ( BGMagenta          "${Esc}[45m" )
    set ( BGCyan             "${Esc}[46m" )
    set ( BGLightGray        "${Esc}[47m" )
    set ( BGDarkGray         "${Esc}[100m" )
    set ( BGLightRed         "${Esc}[101m" )
    set ( BGLightGreen       "${Esc}[102m" )
    set ( BGLightYellow      "${Esc}[103m" )
    set ( BGLightBlue        "${Esc}[104m" )
    set ( BGLightMagenta     "${Esc}[105m" )
    set ( BGLightCyan        "${Esc}[106m" )
    set ( BGWhite            "${Esc}[107m" )

## special colors with 256 color support

    set ( TextAmber          "${Esc}[38;5;208m" )
    set ( BGAmber            "${Esc}[48;5;208m" )

else  ( NOT WIN32 )

    set ( TextReset          "" )

    set ( TextBold           "" )
    set ( TextDim            "" )
    set ( TextItalic         "" )
    set ( TextUnderline      "" )
    set ( TextBlink          "" )
    set ( TextReverse        "" )
    set ( TextHidden         "" )

    set ( TextBoldReset      "" )
    set ( TextDimReset       "" )
    set ( TextItalicReset    "" )
    set ( TextUnderlineReset "" )
    set ( TextBlinkReset     "" )
    set ( TextReverseReset   "" )
    set ( TextHiddenReset    "" )

    set ( TextColorReset     "" )
    set ( TextBlack          "" )
    set ( TextRed            "" )
    set ( TextGreen          "" )
    set ( TextYellow         "" )
    set ( TextBlue           "" )
    set ( TextMagenta        "" )
    set ( TextCyan           "" ) 
    set ( TextLightGray      "" )
    set ( TextDarkGray       "" )
    set ( TextLightRed       "" )
    set ( TextLightGreen     "" )
    set ( TextLightYellow    "" )
    set ( TextLightBlue      "" )
    set ( TextLightMagenta   "" )
    set ( TextLightCyan      "" )
    set ( TextWhite          "" )

    set ( BGColorReset       "" )

    set ( BGBlack            "" )
    set ( BGRed              "" )
    set ( BGGreen            "" )
    set ( BGYellow           "" )
    set ( BGBlue             "" )
    set ( BGMagenta          "" )
    set ( BGCyan             "" )
    set ( BGLightGray        "" )
    set ( BGDarkGray         "" )
    set ( BGLightRed         "" )
    set ( BGLightGreen       "" )
    set ( BGLightYellow      "" )
    set ( BGLightBlue        "" )
    set ( BGLightMagenta     "" )
    set ( BGLightCyan        "" )
    set ( BGWhite            "" )

## special colors with 256 color support

    set ( TextAmber          "" )
    set ( BGAmber            "" )

endif ( NOT WIN32 )

## Help function that prints colored text messages
## inspired by soci colormsg function

function    ( formatText )
    # first arg is name of "return variable"
    list ( GET ARGV 0 RESULT_NAME )
    list ( REMOVE_AT ARGV 0 )

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

    if    ( CMAKE_COLOR_MAKEFILE )
        set ( str "${str}${TextReset}${BGReset}" )
    endif ( CMAKE_COLOR_MAKEFILE )  
    set ( ${RESULT_NAME} ${str} PARENT_SCOPE )
endfunction ( formatText )

function    ( emptyline )
    message ( STATUS "" )
endfunction ( emptyline )

function    ( indent_message LEVEL MSG )

    math ( EXPR NUM_BLANKS "2 * (${LEVEL} - 1)" )
    createBlanks ( TYPE_INTENT ${NUM_BLANKS} )

    message ( STATUS "${TYPE_INTENT}${MSG}" )
endfunction ( indent_message LEVEL MSG )

function    ( heading TEXT )
    emptyline()
    formatText ( H1 "TextUnderline" "${TEXT}" )
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
        formatText ( H2 "${FLAG_FORMAT}" "${FLAG_TEXT}" )
    endif ( VAR STREQUAL "" )
        
    indent_message ( "2" "${TEXT} ${H2}" )
    #emptyline()
endfunction ( heading2 TEXT )

function    ( heading3 TEXT VAR )
    emptyline()
    if    ( VAR STREQUAL "" )
        set ( H2 "" )
    else ( VAR STREQUAL "" )
        if    ( ${VAR} )
            set ( FLAG_TEXT "ENABLED" )
            set ( FLAG_FORMAT "TextGreen" )
        else  ( ${VAR} )
            set ( FLAG_TEXT "DISABLED" )
            set ( FLAG_FORMAT "TextAmber" )
        endif ( ${VAR} )
        formatText ( H3 "${FLAG_FORMAT}" "${FLAG_TEXT}" )
    endif ( VAR STREQUAL "" )
        
    indent_message ( "3" "${TEXT} ${H3}" )
endfunction ( heading3 TEXT )

function    ( found_message TEXT VAR OPTIONAL ADDITIONAL_TEXT )

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
    formatText ( H4 "${FLAG_FORMAT}" "${FLAG_TEXT}" )
        
    indent_message ( "4" "${TEXT}${SCAI_PACKAGE_NAME_BLANKS}${H4} ${ADDITIONAL_TEXT}" )

endfunction ( found_message TEXT )
