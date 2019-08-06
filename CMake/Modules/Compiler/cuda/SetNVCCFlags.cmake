###
 # @file SetNVCCFlags.cmake
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
 # @brief CompilerFlags for LAMA
 # @author Jan Ecker
 # @date 25.04.2013
###

#### Output variables set by this function
# 
#   SCAI_NVCC_FLAGS         : additional flags for NVCC, specific for SCAI projects
#   SCAI_NVCC_FLAGS_RELEASE : additional flags for NVCC, specific for SCAI projects, only Release
#   SCAI_NVCC_FLAGS_DEBUG   : additional flags for NVCC, specific for SCAI projects, only Debug
#
#   The corresponding flags in the variables might be subject of change later by the user

include ( scai_function/checkValue )

#### CUDA specific compiler flags ####

if    ( CUDA_FOUND AND USE_CUDA )
    
    ### choosing the right compute capability
    ### we just start from version 1.3 ( 1.0 - 1.2 is not supported )
    list ( APPEND CC_CHOICES "not-found" "13" "20" "21" "30" "32" "35" "37" "50" "52" "53" "60" "61" "62" "70" "75" "80" )
    set ( CACHE CUDA_COMPUTE_CAPABILITY PROPERTY STRINGS ${CC_CHOICES} )
    checkValue( ${CUDA_COMPUTE_CAPABILITY} "${CC_CHOICES}" )
    set ( CUDA_COMPUTE_CAPABILITY ${CUDA_COMPUTE_CAPABILITY} CACHE STRING "CUDA compute capability (supported up from 13)" )
    
    if    ( CUDA_COMPUTE_CAPABILITY STREQUAL "not-found" )
        set ( CUDA_HAVE_GPU FALSE )
    else  ( CUDA_COMPUTE_CAPABILITY STREQUAL "not-found" )
        set ( CUDA_HAVE_GPU TRUE )
    endif ( CUDA_COMPUTE_CAPABILITY STREQUAL "not-found" )
    mark_as_advanced ( CUDA_COMPUTE_CAPABILITY )

    set ( CUDA_VERBOSE_BUILD OFF )
    set ( CUDA_BUILD_EMULATION OFF )
    
    # unfortunately we can not propagate the host flags to CUDA
    # because this issues to much warning in cuda headers
    # TODO: maybe we can change this with future CUDA releases
    if    ( WIN32 )
        set ( CUDA_PROPAGATE_HOST_FLAGS ON )
    else  ( WIN32 )
        set ( CUDA_PROPAGATE_HOST_FLAGS OFF )
        
        set ( SCAI_NVCC_FLAGS -Xcompiler -fPIC )
        set ( SCAI_NVCC_FLAGS_DEBUG -g -G )
        set ( SCAI_NVCC_FLAGS_RELEASE -O3 -use_fast_math -Xcompiler -ffast-math -Xcompiler -fno-inline )
        # Note: -Xcompiler;-fno-inline is used because of compability issues of CUDA with gcc-4.4

        if    ( CXX_SUPPORTS_C11 )
            if ( CUDA_VERSION_MAJOR LESS 7 )
                message ( FATAL_ERROR "CUDA version ${CUDA_VERSION} does not support -std=c++11, upgrade to CUDA 7.0 or higher" )
            else ()
                list ( APPEND SCAI_NVCC_FLAGS "-std=c++11" )
                if ( CMAKE_CXX_COMPILER_ID MATCHES Clang )
                    list ( APPEND SCAI_NVCC_FLAGS -Xcompiler -stdlib=libc++ )
                endif ( CMAKE_CXX_COMPILER_ID MATCHES Clang )
            endif ()
        endif ( CXX_SUPPORTS_C11 )

        # Intel compiler
        if    ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
            list ( APPEND SCAI_NVCC_FLAGS --compiler-bindir ${CMAKE_CXX_COMPILER}; )  
        endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
        
        # set -march=core02,-mmmx,-msse,-msse2,-msse3,-mssse3,-msse4a flags here
        if    ( MARCH_NATIVE_SUPPORT )
            list ( APPEND SCAI_NVCC_FLAGS_RELEASE -Xcompiler -march=native )
        endif ( MARCH_NATIVE_SUPPORT )
        
    endif ( WIN32 )
    
    if ( NOT CUDA_cusparse_LIBRARY )
        ### cusparse is usually in same directory as cublas
        get_filename_component( HINT_CUDA_LIBRARY_DIR ${CUDA_cublas_LIBRARY} PATH )
        find_library( CUDA_cusparse_LIBRARY NAMES cusparse
                      HINTS ${HINT_CUDA_LIBRARY_DIR} )
        mark_as_advanced( CUDA_cusparse_LIBRARY )
    endif ( NOT CUDA_cusparse_LIBRARY )
    
endif ( CUDA_FOUND AND USE_CUDA )

if ( CUDA_FOUND  )
    set ( CUDA_NVCC_FLAGS_MINSIZEREL     "${CUDA_NVCC_FLAGS_MINSIZEREL}"     CACHE INTERNAL "" )
    set ( CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO}" CACHE INTERNAL "" )
    set ( CUDA_GENERATED_OUTPUT_DIR      "${CUDA_GENERATED_OUTPUT_DIR}"      CACHE INTERNAL "" )
    set ( CUDA_SDK_ROOT_DIR              "${CUDA_SDK_ROOT_DIR}"              CACHE INTERNAL "" )
endif ( CUDA_FOUND  )
