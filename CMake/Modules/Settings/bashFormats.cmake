###
 # @file CMake/Modules/Settings/bashFormats.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @endlicense
 #
 # @brief Sets all escape sequences for bash formatting
 # @author Lauretta Schubert
 # @date 06.04.2016
###

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