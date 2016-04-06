###
 # @file bashFormats.cmake
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
 # @brief Sets all escape sequences for bash formatting
 # @author Lauretta Schubert
 # @date 06.04.2016
 # @since 2.0.0
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

    set ( TextResetBold      "${Esc}[21m" )
    set ( TextResetDim       "${Esc}[22m" )
    set ( TextResetItalic    "${Esc}[23m")
    set ( TextResetUnderline "${Esc}[24m" )
    set ( TextResetBlink     "${Esc}[25m" )
    set ( TextResetReverse   "${Esc}[27m" )
    set ( TextResetHidden    "${Esc}[28m" )

    set ( TextColorReset    "${Esc}[m" )
    set ( TextBlack         "${Esc}[30m" )
    set ( TextRed           "${Esc}[31m" )
    set ( TextGreen         "${Esc}[32m" )
    set ( TextYellow        "${Esc}[33m" )
    set ( TextBlue          "${Esc}[34m" )
    set ( TextMagenta       "${Esc}[35m" )
    set ( TextCyan          "${Esc}[36m" )  
    set ( TextLightGray     "${Esc}[37m" )
    set ( TextDarkGray      "${Esc}[90m" )
    set ( TextLightRed      "${Esc}[91m" )
    set ( TextLightGreen    "${Esc}[92m" )
    set ( TextLightYellow   "${Esc}[93m" )
    set ( TextLightBlue     "${Esc}[94m" )
    set ( TextLightMagenta  "${Esc}[95m" )
    set ( TextLightCyan     "${Esc}[96m" )
    set ( TextWhite         "${Esc}[97m" )

    set ( BGColorReset      "${Esc}[49m" )

    set ( BGBlack           "${Esc}[40m" )
    set ( BGRed             "${Esc}[41m" )
    set ( BGGreen           "${Esc}[42m" )
    set ( BGYellow          "${Esc}[43m" )
    set ( BGBlue            "${Esc}[44m" )
    set ( BGMagenta         "${Esc}[45m" )
    set ( BGCyan            "${Esc}[46m" )
    set ( BGLightGray       "${Esc}[47m" )
    set ( BGDarkGray        "${Esc}[100m" )
    set ( BGLightRed        "${Esc}[101m" )
    set ( BGLightGreen      "${Esc}[102m" )
    set ( BGLightYellow     "${Esc}[103m" )
    set ( BGLightBlue       "${Esc}[104m" )
    set ( BGLightMagenta    "${Esc}[105m" )
    set ( BGLightCyan       "${Esc}[106m" )
    set ( BGWhite           "${Esc}[107m" )

## special colors with 256 color support

    set ( TextAmber         "${Esc}[38;5;208m" )
    set ( BGAmber           "${Esc}[48;5;208m" )

else  ( NOT WIN32 )

    # Todo: fill with empty strings

endif ( NOT WIN32 )