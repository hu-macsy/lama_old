 ###
 # @file FindFortranBLAS.cmake
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
 # @brief Find Fortran BLAS
 # @author
 # @date 25.04.2013
 # $Id$
###

#find_package( BLAS )
find_package( LAPACK ${LAMA_FIND_PACKAGE_FLAGS} )

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
