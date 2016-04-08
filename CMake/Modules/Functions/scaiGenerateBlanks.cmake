###
 # @file Functions/scaiGenerateBlanks.cmake
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
 # @brief Simple internal helper function that generates a blank string that fits the size of an given STRING to LENGTH
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

function    ( createBlanks OUTPUT LENGTH )
	set ( MESSAGE_BLANKS "")
    foreach    ( SCAI_I RANGE ${LENGTH} )
        set ( MESSAGE_BLANKS "${MESSAGE_BLANKS} " )
    endforeach ( SCAI_I RANGE ${SCAI_MESSAGE_BLANK_LENGTH} )

    set ( ${OUTPUT} ${MESSAGE_BLANKS} PARENT_SCOPE )
endfunction ( createBlanks OUTPUT LENGTH )

function    ( scai_generate_blanks OUTPUT STRING LENGTH )
    string ( LENGTH "${STRING}" SCAI_STRING_LENGTH )
    # -1 for correct looping from 0 to LENGTH
    math ( EXPR SCAI_MESSAGE_BLANK_LENGTH ${LENGTH}-${SCAI_STRING_LENGTH} )
    
    createBlanks ( MESSAGE_BLANKS ${SCAI_MESSAGE_BLANK_LENGTH} )
    
    set ( ${OUTPUT} ${MESSAGE_BLANKS} PARENT_SCOPE )
endfunction ( scai_generate_blanks )
