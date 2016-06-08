###
 # @file Functions/setAndCheckCache.cmake
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
 # @brief CMake functions and macros
 # @author Jan Ecker
 # @date 25.04.2013
###

# Function for setting USE_{PACKAGE_NAME} variables depending on {PACKAGE_NAME}_FOUND.
# Also sets cache Variables
function    ( setAndCheckCache PACKAGE_NAME )
	
    # if optional parameter is set, use this one as package name
    if    ( DEFINED ARGV1 )
        set ( CACHE_NAME ${ARGV1} )
    else  ( DEFINED ARGV1 )
    	set ( CACHE_NAME ${PACKAGE_NAME} )
    endif ( DEFINED ARGV1 )

    # Create variable names with USE_XXX and FOUND_XXX
    set ( CACHE_VARIABLE_NAME USE_${CACHE_NAME} )
    set ( FOUND_VARIABLE_NAME ${PACKAGE_NAME}_FOUND )

    # Check if cache variable is already set
    if    ( DEFINED ${CACHE_VARIABLE_NAME} )
        # do nothing
    # if cache variable is NOT set
    else ( DEFINED ${CACHE_VARIABLE_NAME} )
        # Check if package was found
        if    ( ${FOUND_VARIABLE_NAME} )
            set ( USE_PACKAGE TRUE )
        else  ( ${FOUND_VARIABLE_NAME} )
            set ( USE_PACKAGE FALSE )
        endif ( ${FOUND_VARIABLE_NAME} )
        
        # Set cache variable
        set ( ${CACHE_VARIABLE_NAME} ${USE_PACKAGE} CACHE BOOL "" )
    endif ( DEFINED ${CACHE_VARIABLE_NAME} )
endfunction ( setAndCheckCache )
