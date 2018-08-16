###
 # @file FindFortranBLAS.cmake
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
 # @brief Find Fortran BLAS
 # @author 
 # @date 25.04.2013
###

#find_package( BLAS )
find_package( LAPACK ${LAMA_FIND_PACKAGE_FLAGS} )

if    ( APPLE )
	if      ( BLAS_blas_LIBRARY AND LAPACK_lapack_LIBRARY )
	    # everything ok, nothing to be done
	elseif      ( EXISTS ${BLAS_LIBRARIES}/Accelerate )
	    set ( BLAS_blas_LIBRARY ${BLAS_LIBRARIES}/Accelerate )
	else    ()
	    set ( BLAS_FOUND FALSE )
	    set ( LAPACK_FOUND FALSE )
	endif   ()
endif ( APPLE )

## hide flags (we do not use) from the default CMake screen 

set ( BLAS_Accelerate_LIBRARY "${BLAS_Accelerate_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_acml_LIBRARY "${BLAS_acml_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_acml_mp_LIBRARY "${BLAS_acml_mp_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_cblas_LIBRARY "${BLAS_cblas_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_complib.sgimath_LIBRARY "${BLAS_complib.sgimath_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_cxml_LIBRARY "${BLAS_cxml_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_dxml_LIBRARY "${BLAS_dxml_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_essl_LIBRARY "${BLAS_essl_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_scsl_LIBRARY "${BLAS_scsl_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_sgemm_LIBRARY "${BLAS_sgemm_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_sunperf_LIBRARY "${BLAS_sunperf_LIBRARY}" CACHE INTERNAL "" )
set ( BLAS_vecLib_LIBRARY "${BLAS_vecLib_LIBRARY}" CACHE INTERNAL "" )

set ( LAPACK_Accelerate_LIBRARY "${LAPACK_Accelerate_LIBRARY}" CACHE INTERNAL "" )
set ( LAPACK_acml_LIBRARY "${LAPACK_acml_LIBRARY}" CACHE INTERNAL "" )
set ( LAPACK_acml_mp_LIBRARY "${LAPACK_acml_mp_LIBRARY}" CACHE INTERNAL "" )
set ( LAPACK_vecLib_LIBRARY "${LAPACK_vecLib_LIBRARY}" CACHE INTERNAL "" )
