###
 # @file FindLAMA_BLAS.cmake
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
 # @brief CMake functions and macros
 # @author
 # @date 25.04.2013
 # $Id$
###

# - Try to find LAMA_BLAS
#   Once done this will define
#   LAMA_BLAS_FOUND - System has LAMA_BLAS
#   LAMA_BLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

# try to find one of this blas libraries in this order: MKL, ACML, GOTOBLAS, FortranBLAS
include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindMKL.cmake )
include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindACML.cmake )
include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindGOTOBLAS.cmake )
include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindFortranBLAS.cmake )

if ( MKL_FOUND )
    set ( LAMA_BLAS_AUTO "MKL")
elseif ( ACML_FOUND )
    set ( LAMA_BLAS_AUTO "ACML")
elseif ( GOTOBLAS_FOUND )
    set ( LAMA_BLAS_AUTO "GOTOBLAS" )
elseif ( BLAS_FOUND )
    set ( LAMA_BLAS_AUTO "BLAS")
endif ()

LIST ( APPEND LIBRARY_CHOICES "auto" "MKL" "ACML" "GOTOBLAS" "BLAS" )
if ( NOT DEFINED LAMA_BLAS_LIBRARY )
    set ( LAMA_BLAS_LIBRARY "auto" CACHE STRING "Choose the used BLAS Library: ${LIBRARY_CHOICES}" )
else ( NOT DEFINED LAMA_BLAS_LIBRARY )
    set ( LAMA_BLAS_LIBRARY ${LAMA_BLAS_LIBRARY} CACHE STRING "Choose the used BLAS Library: ${LIBRARY_CHOICES}" )
endif ( NOT DEFINED LAMA_BLAS_LIBRARY )
set ( CACHE LAMA_BLAS_LIBRARY PROPERTY STRINGS ${LIBRARY_CHOICES} )
checkValue( ${LAMA_BLAS_LIBRARY} "${LIBRARY_CHOICES}" )

# load selected or auto BLAS Library and set blas style (default: LAMA_FORTRAN_BLAS_STYLE_LOWERCASE)
if ( LAMA_BLAS_LIBRARY STREQUAL "MKL" OR ( LAMA_BLAS_LIBRARY STREQUAL "auto" AND MKL_FOUND ))
    if ( MKL_FOUND )
        set ( LAMA_BLAS_NAME "MKL" )
    	set ( LAMA_BLAS_LIBRARIES ${MKL_LIBRARIES} )
    	set ( LAMA_PBLAS_LIBRARIES ${MKL_PLIBRARIES} )
    	# default: LAMA_FORTRAN_BLAS_STYLE_LOWERCASE
    else ()
        lama_status_message ( ERROR "MKL BLAS library selected, but not found!" )
    endif ( MKL_FOUND )
elseif ( LAMA_BLAS_LIBRARY STREQUAL "ACML" OR ( LAMA_BLAS_LIBRARY STREQUAL "auto" AND ACML_FOUND ))
    if ( ACML_FOUND )
        set ( LAMA_BLAS_NAME "ACML" )
    	set ( LAMA_BLAS_LIBRARIES ${ACML_LIBRARIES} )
    	if ( WIN32 )
            add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UPCASE ) # not tested yet
        else ( WIN32 )
            add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
        endif ( WIN32 )
    else ()
        lama_status_message ( ERROR "ACML BLAS library selected, but not found!" )
    endif ( ACML_FOUND )
elseif ( LAMA_BLAS_LIBRARY STREQUAL "GOTOBLAS" OR ( LAMA_BLAS_LIBRARY STREQUAL "auto" AND GOTOBLAS_FOUND ))
    if ( GOTOBLAS_FOUND )
        set ( LAMA_BLAS_NAME "GOTOBLAS" )
       	set ( LAMA_BLAS_LIBRARIES ${GOTOBLAS_LIBRARIES} )
       	add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
   	else ()
   	    lama_status_message ( ERROR "GOTOBLAS library selected, but not found!" )
   	endif ( GOTOBLAS_FOUND )
#BLAS AND LAPACK found for FortranBLAS
elseif ( LAMA_BLAS_LIBRARY STREQUAL "BLAS" OR ( LAMA_BLAS_LIBRARY STREQUAL "auto" AND BLAS_FOUND ))
    if ( BLAS_FOUND )
        set ( LAMA_BLAS_NAME "BLAS " )
        if ( LAPACK_LIBRARIES )
            set ( LAMA_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
        elseif ( LAPACK_lapack_LIBRARY AND BLAS_blas_LIBRARY )
            list ( APPEND LAMA_BLAS_LIBRARIES ${LAPACK_lapack_LIBRARY} )
        endif ()
        add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
    else ()
        lama_status_message ( ERROR "BLAS library selected, but not found!" )
    endif ( BLAS_FOUND )
endif ()

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set LAMA_BLAS_FOUND to TRUE if all listed variables are TRUE
find_package_handle_standard_args ( LAMA_BLAS DEFAULT_MSG LAMA_BLAS_LIBRARIES)

mark_as_advanced ( LAMA_BLAS_LIBRARIES )
mark_as_advanced ( LAMA_PBLAS_LIBRARIES )
