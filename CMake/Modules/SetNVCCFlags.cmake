###
 # @file SetNVCCFlags.cmake
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
 # @brief CompilerFlags for LAMA
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

#### Output variables set by this function
# 
#   SCAI_NVCC_FLAGS         : additional flags for NVCC, specific for SCAI projects
#   SCAI_NVCC_FLAGS_RELEASE : additional flags for NVCC, specific for SCAI projects, only Release
#   SCAI_NVCC_FLAGS_DEBUG   : additional flags for NVCC, specific for SCAI projects, only Debug
#
#   The corresponding flags in the variables might be subject of change later by the user

include ( Functions/checkValue )

#### CUDA specific compiler flags ####

if    ( CUDA_FOUND AND USE_CUDA )
    
    ### choosing the right compute capability
    ### we just start from version 1.3 ( 1.0 - 1.2 is not supported )
    list ( APPEND CC_CHOICES "not-found" "13" "20" "21" "30" "32" "35" "37" "50" "52" )
	set ( CACHE CUDA_COMPUTE_CAPABILITY PROPERTY STRINGS ${CC_CHOICES} )
    checkValue( ${CUDA_COMPUTE_CAPABILITY} "${CC_CHOICES}" )
    
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
            if ( CUDA_VERSION STRLESS "7.0" )
                message ( FATAL_ERROR "CUDA version ${CUDA_VERSION} does not support -std=c++11, please call cmake with -DCXX_SUPPORTS_C11=0" )
            else ()
                list ( APPEND SCAI_NVCC_FLAGS "-std=c++11" )
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

    ### Check for cuSPASE library, Version 2 (since CUDA 5.0)

    if ( CUDA_VERSION_MAJOR MATCHES "5" )

        message( STATUS "Check for cuSPARSE V2 include file in ${CUDA_INCLUDE_DIRS}" )

        set ( CUSPARSE_V2 false )

        foreach( dir "${CUDA_INCLUDE_DIRS}" )
            if ( EXISTS "${dir}/cusparse_v2.h" )
                set ( CUSPARSE_V2 true )
            endif ( EXISTS "${dir}/cusparse_v2.h" )
        endforeach( dir "${CUDA_INCLUDE_DIRS}" )

        if ( CUSPARSE_V2 )
            message( STATUS "cuSPARSE Version 2 is supported and will be used" )
        else( CUSPARSE_V2 )
            message( STATUS "cuSPARSE Version 2 not supported" )
        endif( CUSPARSE_V2 )

    endif ( CUDA_VERSION_MAJOR MATCHES "5" )
    
endif ( CUDA_FOUND AND USE_CUDA )

if ( CUDA_FOUND  )
    set ( CUDA_NVCC_FLAGS_MINSIZEREL "${CUDA_NVCC_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
    set ( CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO}" CACHE INTERNAL "" )
    set ( CUDA_GENERATED_OUTPUT_DIR "${CUDA_GENERATED_OUTPUT_DIR}" CACHE INTERNAL "" )
    set ( CUDA_SDK_ROOT_DIR "$CUDA_SDK_ROOT_DIR" CACHE INTERNAL "" )
endif ( CUDA_FOUND  )
