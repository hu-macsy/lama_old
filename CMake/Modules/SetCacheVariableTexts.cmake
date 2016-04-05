###
 # @file SetCacheVariableTexts.cmake
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
 # @brief Sets all CacheVariables as previous defined in CMake configuration (default, by User or passed through CMake) with cache texts
 # @author Lauretta Schubert
 # @date 01.04.2016
 # @since 2.0.0
###

## need FORCE directive to force overriding cache text when passing arguments, otherwise it only says: 'initial cache'

# build_

if    ( DEFINED BUILD_DOC )
	set ( BUILD_DOC ${BUILD_DOC} CACHE BOOL "Enable / Disable building of doc" FORCE )
endif ( DEFINED BUILD_DOC )

if    ( DEFINED BUILD_EXAMPLES )
	set ( BUILD_EXAMPLES ${BUILD_EXAMPLES} CACHE BOOL "Enable / Disable building of examples" FORCE )
endif ( DEFINED BUILD_EXAMPLES )

if    ( DEFINED BUILD_TEST )
	set ( BUILD_TEST ${BUILD_TEST} CACHE BOOL "Enable / Disable building of tests" FORCE )
endif ( DEFINED BUILD_TEST )

# cmake_

if    ( DEFINED CMAKE_INSTALL_PREFIX )
	set ( CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "Install path prefix, prepended onto install directories." FORCE )
endif ( DEFINED CMAKE_INSTALL_PREFIX )

if    ( DEFINED CMAKE_BUILD_TYPE )
	set ( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: ${CMAKE_BUILD_TYPE_CHOICES}." FORCE )
endif ( DEFINED CMAKE_BUILD_TYPE )

# cuda_

if    ( DEFINED CUDA_COMPUTE_CAPABILITY )
	set ( CUDA_COMPUTE_CAPABILITY ${CUDA_COMPUTE_CAPABILITY} CACHE STRING "CUDA compute capability (supported up from 13)" FORCE )
endif ( DEFINED CUDA_COMPUTE_CAPABILITY )

if    ( DEFINED CUDA_GENERATE_CODE )
	set ( CUDA_GENERATE_CODE ${CUDA_GENERATE_CODE} CACHE STRING "Which GPU architectures to generate code for (each arch/code pair will be passed as --generate-code option to nvcc, separate multiple pairs by ;)" FORCE )
endif ( DEFINED CUDA_GENERATE_CODE )

# scai_

if    ( DEFINED SCAI_ASSERT_LEVEL )
	set ( SCAI_ASSERT_LEVEL ${SCAI_ASSERT_LEVEL} CACHE STRING "Choose level of ASSERT: ${SCAI_ASSERT_CHOICES}" FORCE )
endif ( DEFINED SCAI_ASSERT_LEVEL )

if    ( DEFINED SCAI_BLAS_LIBRARY )
	set ( SCAI_BLAS_LIBRARY ${SCAI_BLAS_LIBRARY} CACHE STRING "Choose the used BLAS Library: ${SCAI_BLAS_LIBRARY_CHOICES}" FORCE )
endif ( DEFINED SCAI_BLAS_LIBRARY )

if    ( DEFINED SCAI_DOC_TYPE )
	set ( SCAI_DOC_TYPE ${SCAI_DOC_TYPE} CACHE STRING "Choose the type of documentation, options are: ${SCAI_DOC_TYPE_CHOICES}." FORCE )
endif ( DEFINED SCAI_DOC_TYPE )

if    ( DEFINED SCAI_LIBRARY_TYPE )
	set ( SCAI_LIBRARY_TYPE ${SCAI_LIBRARY_TYPE} CACHE STRING "Choose the type of linking: ${SCAI_LIBRARY_TYPE_CHOICES}" FORCE )
endif ( DEFINED SCAI_LIBRARY_TYPE )

if    ( DEFINED SCAI_LOGGING_LEVEL )
	set ( SCAI_LOGGING_LEVEL ${SCAI_LOGGING_LEVEL} CACHE STRING "Choose level of compile time logging: ${SCAI_LOGGING_CHOICES}" FORCE )
endif ( DEFINED SCAI_LOGGING_LEVEL )

if    ( DEFINED SCAI_TRACING )
	set ( SCAI_TRACING ${SCAI_TRACING} CACHE BOOL "Enable / Disable tracing of regions for performance analysis" FORCE )
endif ( DEFINED SCAI_TRACING )

# use_

if    ( DEFINED USE_CODE_COVERAGE )
	set ( USE_CODE_COVERAGE ${USE_CODE_COVERAGE} CACHE BOOL "Enable / Disable use of Code Coverage" FORCE )
endif ( DEFINED USE_CODE_COVERAGE )

if    ( DEFINED USE_CUDA )
	set ( USE_CUDA ${USE_CUDA} CACHE BOOL "Enable / Disable use of CUDA" FORCE )
endif ( DEFINED USE_CUDA )

if    ( DEFINED USE_GPI )
	set ( USE_GPI ${USE_GPI} CACHE BOOL "Enable / Disable use of GPI" FORCE )
endif ( DEFINED USE_GPI )

if    ( DEFINED USE_GRAPHPARTITIONING )
	set ( USE_GRAPHPARTITIONING ${USE_GRAPHPARTITIONING} CACHE BOOL "Enable / Disable use of Graphpartitioning" FORCE )
endif ( DEFINED USE_GRAPHPARTITIONING )

if    ( DEFINED USE_MIC )
	set ( USE_MIC ${USE_MIC} CACHE BOOL "Enable / Disable use of Intel Xeon Phi (MIC)" FORCE )
endif ( DEFINED USE_MIC )

if    ( DEFINED USE_MPI )
	set ( USE_MPI ${USE_MPI} CACHE BOOL "Enable / Disable use of MPI" FORCE )
endif ( DEFINED USE_MPI )

if    ( DEFINED USE_OPENMP )
	set ( USE_OPENMP ${USE_OPENMP} CACHE BOOL "Enable / Disable use of OpenMP" FORCE )
endif ( DEFINED USE_OPENMP )
