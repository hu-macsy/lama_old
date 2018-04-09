###
 # @file scai_function/checkValue.cmake
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
    # need to copy the variable: see CMake docu for this reason
    #Note that the parameters to a macro and values such as ARGN are not variables in the usual CMake sense.
    #They are string replacements much like the C preprocessor would do with a macro.
    set ( MACRO_COPY_OPTION_LIST "${OPTIONLIST}" )
    foreach    ( ITEM ${LIST} )
    	list ( FIND MACRO_COPY_OPTION_LIST ${ITEM} INDEX )
	    if    ( ${INDEX} EQUAL -1 )
            listToString ( ", " "${OPTIONLIST}" ALL_OPTIONS )
	        message ( FATAL_ERROR "Value ${ITEM} is no valid choice out of ${ALL_OPTIONS}" )
	    endif ( ${INDEX} EQUAL -1 )
    endforeach ( ITEM ${VALUELIST} )
endfunction ( checkValues LIST OPTIONLIST )
