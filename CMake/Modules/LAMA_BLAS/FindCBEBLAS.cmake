 ###
 # @file FindCBEBLAS.cmake
 #
 # @license
 # Copyright (c) 2013
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
 # @brief Find CBE BLAS
 # @author
 # @date 25.04.2013
 # $Id$
###

# - Try to find CBE BLAS Implementation 
#
# User can define CBEBLAS_LIBRARY_PATH
#
# Once done this will define 
# CBEBLAS_FOUND - System has LAMA_BLAS
# CBEBLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

## Search CBEBLAS library

find_library ( CBEBLAS_LIBRARY blas HINTS ${CBEBLAS_LIBRARY_PATH} )
if ( EXISTS ${CBEBLAS_LIBRARY} )
    set ( CBEBLAS_LIBRARIES ${CBEBLAS_LIBRARY} )
elseif ( NOT DEFINED CBEBLAS_LIBRARY_PATH )
    message ( STATUS "WARNING CBEBLAS not found. Please define CBEBLAS_LIBRARY_PATH." )
else ()
    message ( STATUS "WARNING CBEBLAS not found. CBEBLAS_LIBRARY_PATH=${CBEBLAS_LIBRARY_PATH} directory does not exist." )
endif ()

## Module footer

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set CBEBLAS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args ( CBEBLAS DEFAULT_MSG CBEBLAS_LIBRARIES )

if( CBEBLAS_FOUND )
   add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
endif( CBEBLAS_FOUND )

mark_as_advanced( CBEBLAS_LIBRARIES CBEBLAS_LIBRARY )