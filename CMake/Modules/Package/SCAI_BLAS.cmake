###
 # @file SCAI_BLAS.cmake
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Find SCAI_BLAS: either MKL, BLAS or use INTERNALBLAS
 # @author 
 # @date 25.04.2013
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

scai_pragma_once ()

#   SCAI_BLAS_FOUND            - System has SCAI_BLAS
#   SCAI_BLAS_NAME             - name of choosen BLAS library
#   SCAI_SCAI_BLAS_INCLUDE_DIR - SCAI_BLAS include directory 
#   SCAI_SCAI_BLAS_LIBRARIES   - The libraries needed to use SCAI_BLAS

enable_language ( C )

if ( NOT DEFINED LAST_SCAI_BLAS_LIBRARY )

    scai_build_variable( NAME      SCAI_BLAS_LIBRARY
                         CHOICES   "auto" "MKL" "BLAS" "INTERNALBLAS" 
                         DEFAULT   "auto"
                         DOCSTRING "Choose the used BLAS library" )

else ( NOT DEFINED LAST_SCAI_BLAS_LIBRARY )

    #message ( STATUS "Choosen BLAS_LIBRARY: ${SCAI_BLAS_LIBRARY}" )

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
        message ( STATUS "MKL BLAS library selected, but not found!" )
    endif ( NOT MKL_FOUND )
elseif ( SCAI_BLAS_LIBRARY STREQUAL "BLAS" )
    find_package ( FortranBLAS ${SCAI_FIND_PACKAGE_FLAGS} )
    if ( NOT BLAS_FOUND )
        message( STATUS "BLAS library selected, but not found!" )    
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
        elseif ( LAPACK_LIBRARIES AND NOT APPLE )
            set ( SCAI_SCAI_BLAS_LIBRARIES ${LAPACK_LIBRARIES} )
	elseif ( BLAS_blas_LIBRARY )
	    set ( SCAI_SCAI_BLAS_LIBRARIES ${BLAS_blas_LIBRARY} )
        endif  ( )
    endif ( BLAS_FOUND )

#message ( STATUS "SCAI_SCAI_BLAS_LIBRARIES ${SCAI_SCAI_BLAS_LIBRARIES}" )

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

if ( "${SCAI_BLAS_NAME}" STREQUAL "MKL" )

    scai_summary_external ( NAME      "SCAI BLAS"
                            ENABLED   True
                            FOUND     True
                            VERSION   "MKL ${MKL_VERSION}"
                            INCLUDE   ${SCAI_SCAI_BLAS_INCLUDE_DIR}
                            LIBRARIES ${SCAI_SCAI_BLAS_LIBRARIES} )

elseif ( "${SCAI_BLAS_NAME}" STREQUAL "BLAS" )

    scai_summary_external ( NAME      "SCAI BLAS"
                            ENABLED   True
                            FOUND     True
                            VERSION   "BLAS ${BLAS_VERSION} Lapack ${LAPACK_VERSION}"
                            INCLUDE   ${SCAI_SCAI_BLAS_INCLUDE_DIR}
                            LIBRARIES ${SCAI_SCAI_BLAS_LIBRARIES} )
else ()

    scai_summary_external ( NAME      "SCAI BLAS"
                            ENABLED   True
                            FOUND     True
                            VERSION   "Internal BLAS"
                          )
endif ()

