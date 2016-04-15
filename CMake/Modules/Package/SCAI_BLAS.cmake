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
 # @brief Find SCAI_BLAS: either MKL, BLAS or use INTERNALBLAS
 # @author
 # @date 25.04.2013
 # @since 1.0.0
###

#   SCAI_BLAS_FOUND            - System has SCAI_BLAS
#   SCAI_BLAS_NAME             - name of choosen BLAS library
#   SCAI_SCAI_BLAS_INCLUDE_DIR - SCAI_BLAS include directory 
#   SCAI_SCAI_BLAS_LIBRARIES   - The libraries needed to use SCAI_BLAS

include ( Functions/checkValue )

enable_language ( C )

if ( NOT DEFINED LAST_SCAI_BLAS_LIBRARY )

if ( NOT DEFINED SCAI_BLAS_LIBRARY )
    set ( SCAI_BLAS_LIBRARY "auto" )
endif ( NOT DEFINED SCAI_BLAS_LIBRARY )
set ( CACHE SCAI_BLAS_LIBRARY PROPERTY STRINGS ${SCAI_BLAS_LIBRARY_CHOICES} )
checkValue( ${SCAI_BLAS_LIBRARY} "${SCAI_BLAS_LIBRARY_CHOICES}" )
set ( SCAI_BLAS_LIBRARY ${SCAI_BLAS_LIBRARY} CACHE STRING "Choose the used BLAS Library: ${SCAI_BLAS_LIBRARY_CHOICES}" )

else ( NOT DEFINED LAST_SCAI_BLAS_LIBRARY )

#message ( STATUS "Choosen BLAS_LIBRARY: ${SCAI_BLAS_LIBRARY}" )

#if ( DEFINED LAST_SCAI_BLAS_LIBRARY )
    if ( NOT LAST_SCAI_BLAS_LIBRARY STREQUAL SCAI_BLAS_LIBRARY )
        set ( SCAI_BLAS_FOUND FALSE )
        set ( SCAI_BLAS_NAME "" )
        set ( SCAI_SCAI_BLAS_LIBRARIES )
      	set ( SCAI_PBLAS_LIBRARIES )
        
        set ( MKL_FOUND FALSE )
        set ( ACML_FOUND FALSE )
        set ( GOTOBLAS_FOUND FALSE )
        set ( BLAS_FOUND FALSE )
        
    endif ( NOT LAST_SCAI_BLAS_LIBRARY STREQUAL SCAI_BLAS_LIBRARY )
#endif ( DEFINED LAST_SCAI_BLAS_LIBRARY )
endif ( NOT DEFINED LAST_SCAI_BLAS_LIBRARY )

# try to find one of this blas libraries in this order: MKL, ACML, GOTOBLAS, FortranBLAS
if ( SCAI_BLAS_LIBRARY STREQUAL "auto" )
    find_package ( MKL ${SCAI_FIND_PACKAGE_FLAGS} )
    if ( NOT MKL_FOUND )
        find_package ( FortranBLAS ${SCAI_FIND_PACKAGE_FLAGS} )
        if ( NOT BLAS_FOUND )
        	set ( INTERNALBLAS_FOUND TRUE )
        endif ( NOT BLAS_FOUND )
    endif ( NOT MKL_FOUND ) 
elseif ( SCAI_BLAS_LIBRARY STREQUAL "MKL" )
    find_package ( MKL ${SCAI_FIND_PACKAGE_FLAGS} )
    if ( NOT MKL_FOUND )
        lama_status_message ( ERROR "MKL BLAS library selected, but not found!" )
    endif ( NOT MKL_FOUND )
elseif ( SCAI_BLAS_LIBRARY STREQUAL "BLAS" )
    find_package ( FortranBLAS ${SCAI_FIND_PACKAGE_FLAGS} )
    if ( NOT BLAS_FOUND )
        lama_status_message ( ERROR "BLAS library selected, but not found!" )    
    endif ( NOT BLAS_FOUND )
#None found or INTERNALBLAS choosen
elseif ( SCAI_BLAS_LIBRARY STREQUAL "INTERNALBLAS")
    set ( INTERNALBLAS_FOUND TRUE ) 
else ()
	message( FATAL_ERROR "none of the possible Library choices (${SCAI_BLAS_LIBRARY_CHOICES}) selected" )
endif ()

if ( NOT INTERNALBLAS_FOUND )

    if ( MKL_FOUND )
        set ( SCAI_BLAS_FOUND TRUE )
        set ( SCAI_BLAS_NAME "MKL" )
        set ( SCAI_SCAI_BLAS_INCLUDE_DIR ${MKL_INCLUDE_DIRS} )
     	set ( SCAI_SCAI_BLAS_LIBRARIES ${MKL_LIBRARIES} )
    endif ( MKL_FOUND )
    
    if ( BLAS_FOUND )
        set ( SCAI_BLAS_FOUND TRUE )
        set ( SCAI_BLAS_NAME "BLAS" )

        if     ( LAPACK_lapack_LIBRARY AND BLAS_blas_LIBRARY )
            set ( SCAI_SCAI_BLAS_LIBRARIES ${BLAS_blas_LIBRARY} ${LAPACK_lapack_LIBRARY} )
        elseif ( LAPACK_LIBRARIES )
            set ( SCAI_SCAI_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
        endif  ( )
    endif ( BLAS_FOUND )

else ( NOT INTERNALBLAS_FOUND )

    set ( SCAI_BLAS_FOUND TRUE )
    set ( SCAI_BLAS_NAME "Internal" )

endif ( NOT INTERNALBLAS_FOUND )

include ( VersionCheck/BLAS )

#message ( STATUS "SCAI_BLAS_FOUND ${SCAI_BLAS_FOUND} SCAI_BLAS_NAME ${SCAI_BLAS_NAME}" )
#message ( STATUS "SCAI_SCAI_BLAS_INCLUDE_DIR ${SCAI_SCAI_BLAS_INCLUDE_DIR} SCAI_SCAI_BLAS_LIBRARIES ${SCAI_SCAI_BLAS_LIBRARIES}" )

include ( FindPackageHandleStandardArgs )
# handle the QUIETLY and REQUIRED arguments and set SCAI_BLAS_FOUND to TRUE if all listed variables are TRUE
#find_package_handle_standard_args ( SCAI_BLAS DEFAULT_MSG SCAI_SCAI_BLAS_LIBRARIES)

mark_as_advanced ( SCAI_SCAI_BLAS_LIBRARIES )

set ( LAST_SCAI_BLAS_LIBRARY ${SCAI_BLAS_LIBRARY} CACHE INTERNAL "" )
