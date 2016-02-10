###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @brief lama/cuda/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
 # @since 1.0.0
###

set( CLASSES      # .cpp, .hpp

        CUDASettings
        CUDATexture
)

set( CUDA_CLASSES  # .cu, .hpp

        CUDAUtils
        CUDACOOUtils
        CUDACSRUtils
        CUDADIAUtils
        CUDAELLUtils
        CUDAJDSUtils
        CUSparseCSRUtils
)

set ( HEADERS        # .hpp

        CUDATexVector
)

if    ( CUSPARSE_V2 )
	set ( set CLASSES ${CLASSES} CUSparseCSRUtils )
	set ( set HEADERS ${HEADERS} CUSPARSETrait CUSPARSEWrapper )
endif ( CUSPARSE_V2 )

### Add classes and header files ###

scai_project(

    CLASSES ${CLASSES} 
    CUDA_CLASSES ${CUDA_CLASSES}
    HEADERS ${HEADERS}

    ADD_PARENT_SCOPE
)
