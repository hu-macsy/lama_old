###
 # @file checkWhiteList.cmake
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
 # @endlicense
 #
 # @brief CMake function to check white list for containing entry
 # @author Lauretta Schubert
 # @date 19.08.2015
###

# checks white list for containing entry and add entry to argument list ( pass whitelist as "${LIST}" (doublequotes !!!) )
macro    ( check_whitelist VAR WHITELIST_NAME ARGS_LIST )

    foreach    ( WHITELIST_ITEM ${WHITELIST_NAME} )
        
        string( COMPARE EQUAL "${CACHE_VAR}" "${WHITELIST_ITEM}" result )
        if    ( ${result} )
        	#message ( STATUS "${CACHE_VAR} matches ${WHITELIST_ITEM}" )
        
        	set ( KEEP_FLAG FALSE )

	        get_property ( CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE )
	        if     ( CACHE_VAR_TYPE STREQUAL "INTERNAL" )
	        	# nothing to do: KEEP_FLAG FALSE
	        elseif ( CACHE_VAR_TYPE MATCHES "PATH" )
	        	if    ( ${${CACHE_VAR}} MATCHES "NOTFOUND" )
	        		#message ( STATUS "both not found ${CACHE_VAR} is ${CACHE_VAR_TYPE} with value ${${CACHE_VAR}}" )
	        		# nothing to do: KEEP_FLAG FALSE
	        	else  ()
					#message ( STATUS "both add ${CACHE_VAR} is ${CACHE_VAR_TYPE} with value ${${CACHE_VAR}}" )
					set ( KEEP_FLAG TRUE )
				endif ()
	        else   ( )
	        	set ( KEEP_FLAG TRUE )
	        endif  ( )

	        if    ( KEEP_FLAG )
	        	#message ( STATUS "ADD ${CACHE_VAR} is ${CACHE_VAR_TYPE} with value ### ${${CACHE_VAR}} ###" )
	            set( CACHE_VAR_TYPE :${CACHE_VAR_TYPE} )
	            list ( APPEND ${ARGS_LIST} "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        else  ( KEEP_FLAG )
				#message ( STATUS "SKIP ${CACHE_VAR} is ${CACHE_VAR_TYPE} with value ### ${${CACHE_VAR}} ###" )
	        endif ( KEEP_FLAG )
        
        endif ( ${result} )
        
    endforeach ( WHITELIST_ITEM ${WHITELIST_NAME} )

endmacro ( check_whitelist )