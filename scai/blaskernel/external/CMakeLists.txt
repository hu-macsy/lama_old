###
 # @file CMakeLists.txt
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
 # @brief blaskernel/external/CMakeLists.txt
 # @author Eric Schricker
 # @date 11.11.2015
###

### Add classes and header files ###

set( CLASSES
#        BLASWrapperNEW
)

set( HEADERS
        LAPACKWrapper
        LAPACKeWrapper
        BLASWrapper
)

if ( SCAI_BLAS_NAME STREQUAL "MKL" )
    if ( EXISTS "${MKL_INCLUDE_DIR}/mkl_lapacke.h" ) 
        set ( LAPACKeFound TRUE )
    endif ( EXISTS "${MKL_INCLUDE_DIR}/mkl_lapacke.h" )
endif ( SCAI_BLAS_NAME STREQUAL "MKL" )

if ( NOT INTERNALBLAS_FOUND )
    set ( CLASSES
        ${CLASSES}
        BLAS_BLAS1
        BLAS_BLAS2
        BLAS_BLAS3
#        LAPACK_LAPACK
    )
endif ( NOT INTERNALBLAS_FOUND )

if ( ${LAPACKeFound} )
    message ( STATUS "LAPACKe is supported and will be used" ) 
    set ( CLASSES ${CLASSES} LAPACKe_LAPACK )
else ( ${LAPACKeFound} )
	if ( NOT INTERNALBLAS_FOUND )
	    message ( STATUS "LAPACKe is not supported, will use LAPACK via Fortran interface" ) 
	    set ( CLASSES ${CLASSES} LAPACK_LAPACK )
    endif ( NOT INTERNALBLAS_FOUND )
endif ( ${LAPACKeFound} )

scai_project ( CLASSES ${CLASSES} HEADERS ${HEADERS} ADD_PARENT_SCOPE )
