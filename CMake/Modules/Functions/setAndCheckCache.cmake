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
