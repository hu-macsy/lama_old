###
 # @file ComputeCapabilityCheck.cmake
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
 # @brief Detect CUDA Compute Capability with test programm
 # @author Jan Ecker
 # @date 10.01.2014
###

if    ( CUDA_FOUND AND NOT CUDA_COMPUTE_CAPABILITY )

    try_run ( CUDA_RUN_RESULT_VAR CUDA_COMPILE_RESULT_VAR
              ${CMAKE_BINARY_DIR}/Compiler/cuda
              ${CMAKE_MODULE_PATH}/Compiler/cuda/ComputeCapabilityCheck.cpp
              CMAKE_FLAGS 
              -DINCLUDE_DIRECTORIES:STRING=${CUDA_TOOLKIT_INCLUDE}
              -DLINK_LIBRARIES:STRING=${CUDA_CUDART_LIBRARY}
              COMPILE_OUTPUT_VARIABLE CUDA_COMPILE_OUTPUT_VAR
              RUN_OUTPUT_VARIABLE CUDA_RUN_OUTPUT_VAR
            )
    
    # CUDA_COMPILE_RESULT_VAR is TRUE when compile succeeds
    # CUDA_RUN_RESULT_VAR is zero when a GPU is found
    if    ( CUDA_COMPILE_RESULT_VAR AND NOT CUDA_RUN_RESULT_VAR )
        set ( CUDA_HAVE_GPU TRUE CACHE BOOL "Whether CUDA-capable GPU is present" )
        set ( CUDA_COMPUTE_CAPABILITY ${CUDA_RUN_OUTPUT_VAR} )
        #set ( CUDA_GENERATE_CODE "arch=compute_${CUDA_COMPUTE_CAPABILITY},code=sm_${CUDA_COMPUTE_CAPABILITY}"
        #      CACHE STRING "Which GPU architectures to generate code for (each arch/code pair will be passed as
        #      --generate-code option to nvcc, separate multiple pairs by ;)" )
        #mark_as_advanced ( CUDA_GENERATE_CODE )
    else  ( CUDA_COMPILE_RESULT_VAR AND NOT CUDA_RUN_RESULT_VAR )
        set ( CUDA_HAVE_GPU FALSE CACHE BOOL "Whether CUDA-capable GPU is present" )
        set ( CUDA_COMPUTE_CAPABILITY "not-found" )
    endif ( CUDA_COMPILE_RESULT_VAR AND NOT CUDA_RUN_RESULT_VAR )
        
    mark_as_advanced ( CUDA_HAVE_GPU CUDA_COMPUTE_CAPABILITY )
   
endif ( CUDA_FOUND AND NOT CUDA_COMPUTE_CAPABILITY )
