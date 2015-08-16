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
function ( lama_status_message )
    string ( ASCII 27 _escape )
    # ANSI Display Atributes
    set ( ERROR "1\;31" )
    set ( WARNING "33" )
    set ( INFO "2\;32" )
    set ( HEADLINE "4" )
    
    set ( coloron FALSE )
    set ( str "" )
    foreach    ( arg ${ARGV} )
        if    ( DEFINED ${arg} AND CMAKE_COLOR_MAKEFILE )
            set(str "${str}${_escape}[${${arg}}m")
            set(coloron TRUE)
        else  ( DEFINED ${arg} AND CMAKE_COLOR_MAKEFILE )
            set ( str "${str}${arg}" )
            if    ( coloron )
                set ( str "${str}${_escape}[0m" )
                set ( coloron FALSE )
            endif ( coloron )
            set ( str "${str} " )
        endif ( DEFINED ${arg} AND CMAKE_COLOR_MAKEFILE )
    endforeach ( arg ${ARGV} )
    
    message ( STATUS ${str} )
endfunction ( lama_status_message )
