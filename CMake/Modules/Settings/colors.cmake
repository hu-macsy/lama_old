## define colors

if    ( NOT WIN32 )

    string ( ASCII 27 Esc )

    set ( TextReset        "${Esc}[0m" )
    set ( TextBold         "${Esc}[1m" )
    set ( TextLowIntensity "${Esc}[2m" )
    set ( TextItalic       "${Esc}[3m" )
    set ( TextUnderline    "${Esc}[4m" )

    set ( TextColorReset   "${Esc}[m" )
    set ( TextBlack        "${Esc}[30m" )
    set ( TextRed          "${Esc}[31m" )
    set ( TextGreen        "${Esc}[32m" )
    set ( TextYellow       "${Esc}[33m" )
    set ( TextBlue         "${Esc}[34m" )
    set ( TextMagenta      "${Esc}[35m" )
    set ( TextCyan         "${Esc}[36m" )  
    set ( TextLightGray    "${Esc}[37m" )
    set ( TextDarkGray     "${Esc}[90m" )
    set ( TextLightRed     "${Esc}[91m" )
    set ( TextLightGreen   "${Esc}[92m" )
    set ( TextLightYellow  "${Esc}[93m" )
    set ( TextLightBlue    "${Esc}[94m" )
    set ( TextLightMagenta "${Esc}[95m" )
    set ( TextLightCyan    "${Esc}[96m" )
    set ( TextWhite        "${Esc}[97m" )

    set ( BGColorReset     "${Esc}[49m" )

    set ( BGBlack          "${Esc}[40m" )
    set ( BGRed            "${Esc}[41m" )
    set ( BGGreen          "${Esc}[42m" )
    set ( BGYellow         "${Esc}[43m" )
    set ( BGBlue           "${Esc}[44m" )
    set ( BGMagenta        "${Esc}[45m" )
    set ( BGCyan           "${Esc}[46m" )
    set ( BGLightGray      "${Esc}[47m" )
    set ( BGDarkGray       "${Esc}[100m" )
    set ( BGLightRed       "${Esc}[101m" )
    set ( BGLightGreen     "${Esc}[102m" )
    set ( BGLightYellow    "${Esc}[103m" )
    set ( BGLightBlue      "${Esc}[104m" )
    set ( BGLightMagenta   "${Esc}[105m" )
    set ( BGLightCyan      "${Esc}[106m" )
    set ( BGWhite          "${Esc}[107m" )

## special colors with 256 color support

    set ( TextAmber         "${Esc}[38;5;208m" )
    set ( BGAmber           "${Esc}[48;5;208m" )

else  ( NOT WIN32 )

    # Todo: fill with empty strings

endif ( NOT WIN32 )