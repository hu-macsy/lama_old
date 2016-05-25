###
 # @file Functions/scaiGenerateBlanks.cmake
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
 # @brief Simple internal helper function that generates a blank string that fits the size of an given STRING to LENGTH
 # @author Jan Ecker
 # @date 25.04.2013
###

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
