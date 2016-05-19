###
 # @file Functions/checkValue.cmake
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
 # @brief CMake functions and macros
 # @author Jan Ecker
 # @date 25.04.2013
###

#checks whether the given value is in the value list ( pass list as "${LIST}" (doublequotes !!!) )
function ( checkValue SINGLEVALUE VALUELIST )
    set ( BOOLVALUE FALSE )
    foreach    ( ITEM ${VALUELIST} )
        if    ( ${SINGLEVALUE} MATCHES ${ITEM} )
            set ( BOOLVALUE TRUE )
        endif ( ${SINGLEVALUE} MATCHES ${ITEM} )
    endforeach ( ITEM ${VALUELIST} )
    if    ( NOT BOOLVALUE )
        message ( FATAL_ERROR "Selected Value ${SINGLEVALUE} is no valid choice out of ${VALUELIST}" )
    endif ( NOT BOOLVALUE )
endfunction ( checkValue SINGLEVALUE VALUELIST )

function ( checkValues LIST OPTIONLIST )
    foreach    ( ITEM ${LIST} )
    	list ( FIND "${OPTIONLIST}" ${ITEM} BOOLVALUE )
	    if    ( NOT BOOLVALUE )
	        message ( FATAL_ERROR "Value ${ITEM} is no valid choice out of ${OPTIONLIST}" )
	    endif ( NOT BOOLVALUE )
    endforeach ( ITEM ${VALUELIST} )
endfunction ( checkValues LIST OPTIONLIST )
