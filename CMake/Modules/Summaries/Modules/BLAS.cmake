###
 # @file CMake/Modules/Summaries/Modules/BLAS.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Summary concerning SCAI_BLAS.
 # @author Lauretta Schubert
 # @date 11.04.2016
###

# BLAS (Lapack)
if    ( SCAI_BLAS_NAME MATCHES "Internal" )
    found_message ( "BLAS" "SCAI_BLAS_FOUND" "REQUIRED" "${SCAI_BLAS_NAME} Version ${SCAI_BLASKERNEL_VERSION}" )
else  ( SCAI_BLAS_NAME MATCHES "Internal" )
    found_message ( "BLAS" "SCAI_BLAS_FOUND" "REQUIRED" "${SCAI_BLAS_NAME} Version ${BLAS_VERSION} with:" )
endif ( SCAI_BLAS_NAME MATCHES "Internal" )

if    ( SCAI_BLAS_NAME MATCHES "BLAS" )
    message ( STATUS "                                 ${BLAS_blas_LIBRARY}" )
    found_message ( "Lapack" "LAPACK_FOUND" "REQUIRED" "with: ${LAPACK_lapack_LIBRARY} " )
else  ( SCAI_BLAS_NAME MATCHES "BLAS" )
    foreach    ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
        message ( STATUS "                                 ${_B_LIB}" )
    endforeach ( _B_LIB ${SCAI_SCAI_BLAS_LIBRARIES} )
endif ( SCAI_BLAS_NAME MATCHES "BLAS" )