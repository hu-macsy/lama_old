###
 # @file SCAI_BLAS.cmake
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
 # @author
 # @date 25.04.2013
 # @since 1.0.0
###

# - Try to find SCAI_BLAS
#   Once done this will define
#   SCAI_BLAS_FOUND         - System has SCAI_BLAS
#   SCAI_BLAS_INCLUDE_DIR   - SCAI_BLAS include directory 
#   SCAI_BLAS_LIBRARIES     - The libraries needed to use SCAI_BLAS

include ( Functions/checkValue )

enable_language ( C )

LIST ( APPEND LIBRARY_CHOICES "auto" "MKL" "BLAS" "INTERNALBLAS" )
if ( NOT DEFINED SCAI_BLAS_LIBRARY )
    set ( SCAI_BLAS_LIBRARY "auto" CACHE STRING "Choose the used BLAS Library: ${LIBRARY_CHOICES}" )
else ( NOT DEFINED SCAI_BLAS_LIBRARY )
    set ( SCAI_BLAS_LIBRARY ${SCAI_BLAS_LIBRARY} CACHE STRING "Choose the used BLAS Library: ${LIBRARY_CHOICES}" )
endif ( NOT DEFINED SCAI_BLAS_LIBRARY )
set ( CACHE SCAI_BLAS_LIBRARY PROPERTY STRINGS ${LIBRARY_CHOICES} )
checkValue( ${SCAI_BLAS_LIBRARY} "${LIBRARY_CHOICES}" )

if ( DEFINED LAST_SCAI_BLAS_LIBRARY )
    if ( NOT LAST_SCAI_BLAS_LIBRARY STREQUAL SCAI_BLAS_LIBRARY )
        set ( SCAI_BLAS_FOUND FALSE )
        set ( SCAI_BLAS_NAME "" )
        set ( SCAI_BLAS_LIBRARIES )
      	set ( SCAI_PBLAS_LIBRARIES )
        
        set ( MKL_FOUND FALSE )
        set ( ACML_FOUND FALSE )
        set ( GOTOBLAS_FOUND FALSE )
        set ( BLAS_FOUND FALSE )
        
    endif ( NOT LAST_SCAI_BLAS_LIBRARY STREQUAL SCAI_BLAS_LIBRARY )
endif ( DEFINED LAST_SCAI_BLAS_LIBRARY )

# try to find one of this blas libraries in this order: MKL, ACML, GOTOBLAS, FortranBLAS
if ( SCAI_BLAS_LIBRARY STREQUAL "auto" )
    include ( ${CMAKE_MODULE_PATH}/SCAI_BLAS/FindMKL.cmake )
    if ( NOT MKL_FOUND )
        include ( ${CMAKE_MODULE_PATH}/SCAI_BLAS/FindFortranBLAS.cmake )
        if ( NOT BLAS_FOUND )
        	set ( INTERNALBLAS_FOUND TRUE )
        endif ( NOT BLAS_FOUND )
    endif ( NOT MKL_FOUND ) 
elseif ( SCAI_BLAS_LIBRARY STREQUAL "MKL" )
    include ( ${CMAKE_MODULE_PATH}/SCAI_BLAS/FindMKL.cmake )
    if ( NOT MKL_FOUND )
        lama_status_message ( ERROR "MKL BLAS library selected, but not found!" )
    endif ( NOT MKL_FOUND )
elseif ( SCAI_BLAS_LIBRARY STREQUAL "BLAS" )
    include ( ${CMAKE_MODULE_PATH}/SCAI_BLAS/FindFortranBLAS.cmake )
    if ( NOT BLAS_FOUND )
        lama_status_message ( ERROR "BLAS library selected, but not found!" )    
    endif ( NOT BLAS_FOUND )
#None found or INTERNALBLAS choosen
elseif ( SCAI_BLAS_LIBRARY STREQUAL "INTERNALBLAS")
    set ( INTERNALBLAS_FOUND TRUE ) 
else ()
	message( FATAL_ERROR "none of the possible Library choices (${LIBRARY_CHOICES}) selected" )
endif ()

if ( NOT INTERNALBLAS_FOUND )

    # load selected or auto choosen BLAS Library and set blas style (default: SCAI_FORTRAN_BLAS_STYLE_LOWERCASE)
    if ( MKL_FOUND )

        set ( SCAI_BLAS_FOUND TRUE )
        set ( SCAI_BLAS_NAME "MKL" )
        set ( SCAI_BLAS_INCLUDE_DIR ${MKL_INCLUDE_DIRS} )
     	set ( SCAI_BLAS_LIBRARIES ${MKL_LIBRARIES} )
      	set ( SCAI_PBLAS_LIBRARIES ${MKL_PLIBRARIES} )
       	# default: SCAI_FORTRAN_BLAS_STYLE_LOWERCASE
       	
    endif ( MKL_FOUND )
    
    if ( BLAS_FOUND )
    
        set ( SCAI_BLAS_FOUND TRUE )
        set ( SCAI_BLAS_NAME "BLAS " )
        set ( SCAI_BLAS_INCLUDE_DIR ${} ) ##TODO
        if ( LAPACK_LIBRARIES )
            set ( SCAI_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
        elseif ( LAPACK_lapack_LIBRARY AND BLAS_blas_LIBRARY )
            list ( APPEND SCAI_BLAS_LIBRARIES ${LAPACK_lapack_LIBRARY} )
        endif ()
        add_definitions ( -DSCAI_FORTRAN_BLAS_STYLE_UNDERSCORE )
    endif ( BLAS_FOUND )

else ( NOT INTERNALBLAS_FOUND )
    set ( SCAI_BLAS_FOUND TRUE )
    set ( SCAI_BLAS_NAME "Internal" )
    set ( SCAI_BLAS_INCLUDE_DIR ${CMAKE_ROOT_DIR}/scai/lama/cblas/include )
endif ( NOT INTERNALBLAS_FOUND )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set SCAI_BLAS_FOUND to TRUE if all listed variables are TRUE
#find_package_handle_standard_args ( SCAI_BLAS DEFAULT_MSG SCAI_BLAS_LIBRARIES)

mark_as_advanced ( SCAI_BLAS_LIBRARIES SCAI_PBLAS_LIBRARIES )

set ( LAST_SCAI_BLAS_LIBRARY ${SCAI_BLAS_LIBRARY} CACHE INTERNAL "" )
