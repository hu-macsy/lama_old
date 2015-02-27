###
 # @file FindLAMA_BLAS.cmake
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

# - Try to find LAMA_BLAS
#   Once done this will define
#   LAMA_BLAS_FOUND - System has LAMA_BLAS
#   LAMA_BLAS_LIBRARIES - The libraries needed to use LAMA_BLAS

LIST ( APPEND LIBRARY_CHOICES "auto" "MKL" "ACML" "GOTOBLAS" "BLAS" "INTERNALBLAS" )
if ( NOT DEFINED LAMA_BLAS_LIBRARY )
    set ( LAMA_BLAS_LIBRARY "auto" CACHE STRING "Choose the used BLAS Library: ${LIBRARY_CHOICES}" )
else ( NOT DEFINED LAMA_BLAS_LIBRARY )
    set ( LAMA_BLAS_LIBRARY ${LAMA_BLAS_LIBRARY} CACHE STRING "Choose the used BLAS Library: ${LIBRARY_CHOICES}" )
endif ( NOT DEFINED LAMA_BLAS_LIBRARY )
set ( CACHE LAMA_BLAS_LIBRARY PROPERTY STRINGS ${LIBRARY_CHOICES} )
checkValue( ${LAMA_BLAS_LIBRARY} "${LIBRARY_CHOICES}" )

if ( DEFINED LAST_LAMA_BLAS_LIBRARY )
    if ( NOT LAST_LAMA_BLAS_LIBRARY STREQUAL LAMA_BLAS_LIBRARY )
        set ( LAMA_BLAS_FOUND FALSE )
        set ( LAMA_BLAS_NAME "" )
        set ( LAMA_BLAS_LIBRARIES )
      	set ( LAMA_PBLAS_LIBRARIES )
        
        set ( MKL_FOUND FALSE )
        set ( ACML_FOUND FALSE )
        set ( GOTOBLAS_FOUND FALSE )
        set ( BLAS_FOUND FALSE )
        
    endif ( NOT LAST_LAMA_BLAS_LIBRARY STREQUAL LAMA_BLAS_LIBRARY )
endif ( DEFINED LAST_LAMA_BLAS_LIBRARY )

# try to find one of this blas libraries in this order: MKL, ACML, GOTOBLAS, FortranBLAS
if ( LAMA_BLAS_LIBRARY STREQUAL "auto" )
    include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindMKL.cmake )
    if ( NOT MKL_FOUND )
        include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindACML.cmake )
        if ( NOT ACML_FOUND )
            include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindGOTOBLAS.cmake )
            if ( NOT GOTOBLAS_FOUND )
                include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindFortranBLAS.cmake )
                if ( NOT BLAS_FOUND )
                    set ( INTERNALBLAS_FOUND TRUE )
                endif ( NOT BLAS_FOUND )
            endif ( NOT GOTOBLAS_FOUND )
        endif (  NOT ACML_FOUND )
    endif ( NOT MKL_FOUND ) 
elseif ( LAMA_BLAS_LIBRARY STREQUAL "MKL" )
    include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindMKL.cmake )
    if ( NOT MKL_FOUND )
        lama_status_message ( ERROR "MKL BLAS library selected, but not found!" )
    endif ( NOT MKL_FOUND )
elseif ( LAMA_BLAS_LIBRARY STREQUAL "ACML" )
    include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindACML.cmake )
    if ( NOT ACML_FOUND )
        lama_status_message ( ERROR "ACML BLAS library selected, but not found!" )
    endif ( NOT ACML_FOUND )
elseif ( LAMA_BLAS_LIBRARY STREQUAL "GOTOBLAS" )
    include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindGOTOBLAS.cmake )
    if ( NOT GOTO_BLAS_FOUND )
        lama_status_message ( ERROR "GOTOBLAS library selected, but not found!" )
    endif ( NOT GOTO_BLAS_FOUND )
elseif ( LAMA_BLAS_LIBRARY STREQUAL "BLAS" )
    include ( ${CMAKE_MODULE_PATH}/LAMA_BLAS/FindFortranBLAS.cmake )
    if ( NOT BLAS_FOUND )
        lama_status_message ( ERROR "BLAS library selected, but not found!" )    
    endif ( NOT BLAS_FOUND )
#None found or INTERNALBLAS choosen
elseif ( LAMA_BLAS_LIBRARY STREQUAL "INTERNALBLAS")
    set ( INTERNALBLAS_FOUND TRUE ) 
else ()
	message( FATAL_ERROR "none of the possible Library choices (${LIBRARY_CHOICES}) selected" )
endif ()

if ( NOT INTERNALBLAS_FOUND )

    # load selected or auto choosen BLAS Library and set blas style (default: LAMA_FORTRAN_BLAS_STYLE_LOWERCASE)
    if ( MKL_FOUND )
        set ( LAMA_BLAS_FOUND TRUE )
        set ( LAMA_BLAS_NAME "MKL" )
     	set ( LAMA_BLAS_LIBRARIES ${MKL_LIBRARIES} )
      	set ( LAMA_PBLAS_LIBRARIES ${MKL_PLIBRARIES} )
       	# default: LAMA_FORTRAN_BLAS_STYLE_LOWERCASE
    endif ( MKL_FOUND )
    
    if ( ACML_FOUND )
        set ( LAMA_BLAS_FOUND TRUE )
        set ( LAMA_BLAS_NAME "ACML" )
        set ( LAMA_BLAS_LIBRARIES ${ACML_LIBRARIES} )
        if ( WIN32 )
            add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UPCASE ) # not tested yet
        else ( WIN32 )
            add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
        endif ( WIN32 )
    endif ( ACML_FOUND )
    
    if ( GOTOBLAS_FOUND )
        set ( LAMA_BLAS_FOUND TRUE )
        set ( LAMA_BLAS_NAME "GOTOBLAS" )
       	set ( LAMA_BLAS_LIBRARIES ${GOTOBLAS_LIBRARIES} )
       	add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
    endif ( GOTOBLAS_FOUND )
    
    if ( BLAS_FOUND )
        set ( LAMA_BLAS_FOUND TRUE )
        set ( LAMA_BLAS_NAME "BLAS " )
        if ( LAPACK_LIBRARIES )
            set ( LAMA_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
        elseif ( LAPACK_lapack_LIBRARY AND BLAS_blas_LIBRARY )
            list ( APPEND LAMA_BLAS_LIBRARIES ${LAPACK_lapack_LIBRARY} )
        endif ()
        add_definitions ( -DLAMA_FORTRAN_BLAS_STYLE_UNDERSCORE )
    endif ( BLAS_FOUND )

else ( NOT INTERNALBLAS_FOUND )
    set ( LAMA_BLAS_FOUND TRUE )
    set ( LAMA_BLAS_NAME "Internal" )
endif ( NOT INTERNALBLAS_FOUND )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set LAMA_BLAS_FOUND to TRUE if all listed variables are TRUE
#find_package_handle_standard_args ( LAMA_BLAS DEFAULT_MSG LAMA_BLAS_LIBRARIES)

mark_as_advanced ( LAMA_BLAS_LIBRARIES )
mark_as_advanced ( LAMA_PBLAS_LIBRARIES )

set ( LAST_LAMA_BLAS_LIBRARY ${LAMA_BLAS_LIBRARY} CACHE INTERNAL "" )
