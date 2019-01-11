###
 # @file Package/FFTW.cmake
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
 # @brief SCAI wrapper for find_package( PNG )
 # @author Thomas Brandes
 # @date 14.05.2017
###

include ( scai_macro/scai_pragma_once )

scai_pragma_once()

#  CMake provides already a module to find the PNG reference library, use it

find_package( FFTW ${SCAI_FIND_PACKAGE_FLAGS} )

if ( FFTW_FOUND )

else ()

    # Try to use the FFTW wrapper of the MKL libraries

    find_package ( MKL ${SCAI_FIND_PACKAGE_FLAGS} )

    # message ( STATUS "FFTW not found, MKL_FOUND=${MKL_FOUND}" )

    if ( MKL_FOUND )

        # message ( STATUS "Try to use FFTW of MKL" )

        find_path( FFTW_INCLUDE_DIR fftw3.h
                   ${MKL_INCLUDE_DIR}/fftw )

        if ( FFTW_INCLUDE_DIR )
            set ( FFTW_FOUND TRUE )
            set ( FFTW_LIBRARIES ${MKL_LIBRARIES} )
        endif ()

    endif ()
    
endif ()

if ( FFTW_FOUND )

    ## ToDo: get FFTW version


    ## set the SCAI variables to get macros for dependencies working correctly

    set ( SCAI_FFTW_INCLUDE_DIR ${FFTW_INCLUDE_DIR} )
    set ( SCAI_FFTW_LIBRARIES ${FFTW_LIBRARIES} )

endif ()

scai_build_variable ( NAME      USE_FFTW
                      BOOL 
                      DEFAULT   ${FFTW_FOUND}
                      DOCSTRING "use of FFTW library (Discrete Fast Fourier Transform)" )

scai_summary_external ( NAME      FFTW
                        ENABLED   ${USE_FFTW}
                        FOUND     ${FFTW_FOUND} 
                        VERSION   "3.3.7"
                        INCLUDE   ${FFTW_INCLUDE_DIR} 
                        LIBRARIES ${FFTW_LIBRARIES} )

