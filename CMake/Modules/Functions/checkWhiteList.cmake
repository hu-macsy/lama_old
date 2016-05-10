###
 # @file checkWhiteList.cmake
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
        
	        get_property ( CACHE_VAR_TYPE CACHE ${CACHE_VAR} PROPERTY TYPE )
	        #message ( STATUS "CACHE_VAR ${CACHE_VAR}" )
	        if     ( CACHE_VAR_TYPE STREQUAL "INTERNAL" )
	            # skip Variable
	        	# message ( STATUS "Skip Variable ${CACHE_VAR}" )
	        #elseif ( CACHE_VAR_TYPE STREQUAL "UNINITIALIZED" )
	        #    set ( CACHE_VAR_TYPE )
	        #    list ( APPEND ${ARGS_LIST} "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        else   ( )
	            set( CACHE_VAR_TYPE :${CACHE_VAR_TYPE} )
	            list ( APPEND ${ARGS_LIST} "-D${CACHE_VAR}${CACHE_VAR_TYPE}=${${CACHE_VAR}}" )
	        endif  ( )
        
        endif ( ${result} )
        
    endforeach ( WHITELIST_ITEM ${WHITELIST_NAME} )

endmacro ( check_whitelist )