###
 # @file CheckHostCompilerCompatibility.cmake
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
 # @brief Check if Host Compiler and CUDA Installation are compatible with test programm
 # @author Lauretta Schubert
 # @date 02.04.2016
 # @since 2.0.0
###

if    ( CUDA_FOUND )

	message ( STATUS "execute: ${CUDA_NVCC_EXECUTABLE} ${CMAKE_MODULE_PATH}/Compiler/cuda/CheckHostCompilerCompatibility.cu -ccbin=${CMAKE_CXX_COMPILER}")

	execute_process ( COMMAND ${CUDA_NVCC_EXECUTABLE} ${CMAKE_MODULE_PATH}/Compiler/cuda/CheckHostCompilerCompatibility.cu -ccbin=${CMAKE_CXX_COMPILER}
					  RESULT_VARIABLE CUDA_CHECK_COMPILE_RESULT_VAR
					  OUTPUT_VARIABLE CUDA_CHECK_COMPILE_OUTPUT_VAR
					  ERROR_VARIABLE  CUDA_CHECK_COMPILE_ERROR_VAR )

	if    ( CUDA_CHECK_COMPILE_RESULT_VAR )
		message ( STATUS "Your CUDA Installation is not supported by your icc version:\n ${CUDA_CHECK_COMPILE_ERROR_VAR}Disable CUDA support." )
		set ( CUDA_FOUND FALSE )
	endif ( CUDA_CHECK_COMPILE_RESULT_VAR )

endif ( CUDA_FOUND )