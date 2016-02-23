###
 # @file checkWhiteList.cmake
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
 # @brief CMake function to check white list for containing entry
 # @author Lauretta Schubert
 # @date 19.08.2015
 # @since 2.0.0
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