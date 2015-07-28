###
 # @file CUDAComputeCapability.cmake
 #
 # @license
 # Copyright (c) 2009-2014
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
 # @brief Version variable defintions for the used compilers
 # @author Jan Ecker
 # @date 10.01.2014
 # @since 1.1.0
 #
 # Based on code by Florian Rathgeber<florian.rathgeber@gmail.com> on
 # http://www.cmake.org/Bug/print_bug_page.php?bug_id=11767
###

if ( CUDA_FOUND )
    try_run ( RUN_RESULT_VAR COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}
        ${CMAKE_MODULE_PATH}/CudaComputeCapability.c
        CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${CUDA_TOOLKIT_INCLUDE}
        -DLINK_LIBRARIES:STRING=${CUDA_CUDART_LIBRARY}
        COMPILE_OUTPUT_VARIABLE COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE RUN_OUTPUT_VAR)
        
    # COMPILE_RESULT_VAR is TRUE when compile succeeds
    # RUN_RESULT_VAR is zero when a GPU is found
    if ( COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR )
        set ( CUDA_HAVE_GPU TRUE CACHE BOOL "Whether CUDA-capable GPU is present" )
        set ( CUDA_COMPUTE_CAPABILITY ${RUN_OUTPUT_VAR} CACHE STRING "CUDA compute capability (supported up from 13)")
        set ( CUDA_GENERATE_CODE "arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY}" CACHE STRING "Which GPU architectures to generate code for (each arch/code pair will be passed as --generate-code option to nvcc, separate multiple pairs by ;)" )
        mark_as_advanced ( CUDA_COMPUTE_CAPABILITY CUDA_GENERATE_CODE )
    else ( COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR )
        set ( CUDA_HAVE_GPU FALSE CACHE BOOL "Whether CUDA-capable GPU is present")
        set ( CUDA_COMPUTE_CAPABILITY "not-found" CACHE STRING "CUDA compute capability (supported up from 13)" )
    endif ( COMPILE_RESULT_VAR AND NOT RUN_RESULT_VAR )
    mark_as_advanced( CUDA_HAVE_GPU )
    
    if ( CUDA_COMPUTE_CAPABILITY MATCHES "35" )
        add_definitions ( -DCUDA_35 )
    endif ( CUDA_COMPUTE_CAPABILITY MATCHES "35" )
endif ( CUDA_FOUND )