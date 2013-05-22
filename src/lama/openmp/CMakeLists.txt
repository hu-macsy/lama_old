###
 # @file CMakeLists.txt
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
 # @brief lama/openmp/CMakeLists.txt
 # @author Jan Ecker
 # @date 02.04.2013
 # @since 1.0.0
###

### Add classes and header files ###
set ( CLASSES
        BLASHelper
        OpenMPBLAS1
        OpenMPBLAS2
        OpenMPBLAS3
        OpenMPCOOUtils
        OpenMPCSRUtils
        OpenMPDenseUtils
        OpenMPDIAUtils
        OpenMPELLUtils
        OpenMPJDSUtils
        OpenMPUtils
    )

# OpenMPLAPACKe requires LAPACKe interface, is in MKL

set ( LAPACKEeFound FALSE )

if ( LAMA_BLAS_NAME STREQUAL "MKL" )
    if ( EXISTS "${MKL_INCLUDE_DIR}/mkl_lapacke.h" ) 
        set ( LAPACKeFound TRUE )
    endif()
endif()

if ( ${LAPACKeFound} )
    message( STATUS "LAPACKe is supported and will be used" ) 
    set ( CLASSES
            ${CLASSES}
            OpenMPLAPACKe
    )
else()
    message( STATUS "LAPACKe is not supported, will use LAPACK via Fortran interface" ) 
    set ( CLASSES
            ${CLASSES}
            OpenMPLAPACK
    )
endif()
    
lama_classes ( ${CLASSES} )

### add sources and headers to parent scope ###
lama_add ()
