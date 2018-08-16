###
 # @file CheckHostCompilerCompatibility.cmake
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
 # @brief Check if Host Compiler and CUDA Installation are compatible with test programm
 # @author Lauretta Schubert
 # @date 02.04.2016
###

if ( CUDA_FOUND )

    if ( SCAI_CMAKE_VERBOSE )
        message ( STATUS "Check compatibility of ${CUDA_NVCC_EXECUTABLE} with ${CUDA_HOST_COMPILER}" )
    endif ()

    execute_process ( COMMAND ${CUDA_NVCC_EXECUTABLE} ${CMAKE_MODULE_PATH}/Compiler/cuda/CheckHostCompilerCompatibility.cu -ccbin=${CUDA_HOST_COMPILER}
                      RESULT_VARIABLE CUDA_CHECK_COMPILE_RESULT_VAR
                      OUTPUT_VARIABLE CUDA_CHECK_COMPILE_OUTPUT_VAR
                      ERROR_VARIABLE  CUDA_CHECK_COMPILE_ERROR_VAR )

    if ( SCAI_CMAKE_VERBOSE )
        message ( STATUS "Result = ${CUDA_CHECK_COMPILE_RESULT_VAR}" )
    endif ()

    if ( CUDA_CHECK_COMPILE_RESULT_VAR )

        if ( CMAKE_CXX_COMPILER_ID MATCHES GNU )
            string ( REGEX MATCH "gcc ([0-9]+\\.[0-9]) and up" VERSION_NUMBER ${CUDA_CHECK_COMPILE_ERROR_VAR} )
            message ( STATUS "Current CUDA Installation does not work with ${VERSION_NUMBER} - you have GCC ${CXX_COMPILER_VERSION}. Disable CUDA support." )
        endif ( CMAKE_CXX_COMPILER_ID MATCHES GNU )

        if ( CMAKE_CXX_COMPILER_ID MATCHES Clang )
            string ( REGEX MATCH "\\(\\'([0-9]+)\\'\\) of the host compiler" VERSION_NUMBER ${CUDA_CHECK_COMPILE_ERROR_VAR} )
            message ( STATUS "Current CUDA Installation does not work with CLANG ${VERSION_NUMBER}. Disable CUDA support." )
        endif ( CMAKE_CXX_COMPILER_ID MATCHES Clang )

        if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
            string ( REGEX MATCH "([0-9]+\\.[0-9])" VERSION_NUMBER ${CUDA_CHECK_COMPILE_ERROR_VAR} )
            message ( STATUS "Current CUDA Installation only works with ICC ${VERSION_NUMBER} - you have ICC ${CXX_COMPILER_VERSION}. Disable CUDA support." )
        endif ()

        set ( CUDA_FOUND FALSE )

    endif ()

endif ()
