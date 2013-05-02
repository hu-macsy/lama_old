###
 # @file CMakeLists.txt
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
 # @brief lama/cuda/CMakeLists.txt
 # @author Jan Ecker
 # @date 03.04.2013
###

### set this dir as lama source dir ###
lama_set_source_dir ()

### Add classes and header files ###
set ( CLASSES
        CUDABLAS1
        CUDABLAS2
        CUDABLAS3
        CUDAContext
        CUDAContextManager
        CUDAError
        CUDAException
        CUDAHostContext
        CUDAHostContextManager
        CUDALAPACK
        CUDAStreamSyncToken
        CUDATracerSyncToken
    )

set ( HEADERS
        CUDATracerHelper
    )
    
lama_classes ( ${CLASSES} )
lama_headers ( ${HEADERS} )

### Add CUDA Sources ###
set ( CUDA_SOURCES
        blas1.cu
        blas3.cu
        emptyKernel.cu
        CUDAUtils.cu
        CUDACSRUtils.cu
        CUDAELLUtils.cu
        CUDAJDSUtils.cu
        CUDADIAUtils.cu
        CUDACOOUtils.cu
    )
    
### set include directories ### 
include_directories ( . ${CUDA_INCLUDE_DIRS} )

### add CUDA library ###

if ( NOT WIN32 )

    # Linux requires explicit CUDA compilation of *.cu files

    cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
    add_library( lamacuda SHARED ${CXX_SOURCES} ${CUDA_FILES} )

else( NOT WIN32 )

    # Windows suports different compilation for *.cpp and *.cu files

    cuda_add_library( lamacuda SHARED ${CXX_SOURCES} ${CUDA_SOURCES} )

endif( NOT WIN32 )


### target link libraries ###
target_link_libraries( lamacuda
                        ${CUDA_cusparse_LIBRARY}
                        ${CUDA_cublas_LIBRARY} 
                        ${CUDA_CUDART_LIBRARY}
                        ${CUDA_CUDA_LIBRARY}
                     )

### install library
install ( TARGETS lamacuda DESTINATION lib )
