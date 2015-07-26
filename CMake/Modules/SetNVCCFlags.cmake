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

#### CUDA specific compiler flags ####

if ( CUDA_FOUND AND LAMA_USE_CUDA )
    
    ### choosing the right compute capability
    ### we just start from version 1.3 ( 1.0 - 1.2 is not supported )
    LIST ( APPEND CC_CHOICES "not-found" "13" "20" "21" "30" "35" )
    #set ( CUDA_COMPUTE_CAPABILITY "13" CACHE STRING "CUDA compute capability (supported up from 13)" )
	set ( CACHE CUDA_COMPUTE_CAPABILITY PROPERTY STRINGS ${CC_CHOICES} )
    checkValue( ${CUDA_COMPUTE_CAPABILITY} "${CC_CHOICES}" )
    if ( CUDA_COMPUTE_CAPABILITY STREQUAL "not-found" )
        set ( CUDA_HAVE_GPU FALSE )
    else ( CUDA_COMPUTE_CAPABILITY STREQUAL "not-found" )
        set ( CUDA_HAVE_GPU TRUE )
    endif ( CUDA_COMPUTE_CAPABILITY STREQUAL "not-found" )
	mark_as_advanced ( CUDA_COMPUTE_CAPABILITY )

    set ( CUDA_VERBOSE_BUILD OFF )
    set ( CUDA_BUILD_EMULATION OFF )
    
    # unfortunately we can not propagate the host flags to CUDA
    # because this issues to much warning in cuda headers
    # TODO: maybe we can change this with future CUDA releases
    if ( WIN32 )
        set ( CUDA_PROPAGATE_HOST_FLAGS ON )
    else ( WIN32 )
        set ( CUDA_PROPAGATE_HOST_FLAGS OFF )
        
        set ( ADDITIONAL_NVCC_FLAGS -Xcompiler -fPIC )
        set ( ADDITIONAL_NVCC_RELEASE_FLAGS -O3 -use_fast_math -Xcompiler -ffast-math -Xcompiler -fno-inline )
        
        if ( CXX_SUPPORTS_C11 )
            list ( APPEND ADDITIONAL_NVCC_FLAGS "-std=c++11" )
        endif ( CXX_SUPPORTS_C11 )

        # Intel compiler
        if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
            list ( APPEND ADDITIONAL_NVCC_FLAGS "--compiler-bindir ${CMAKE_CXX_COMPILER};" )  
        endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
        
        #-Xcompiler;-fno-inline is used because of compability issues of CUDA with gcc-4.4
        if ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
      	    list ( APPEND ADDITIONAL_NVCC_FLAGS -g -G )
        endif ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
        
        # set -march=core02,-mmmx,-msse,-msse2,-msse3,-mssse3,-msse4a flags here
        if ( MARCH_NATIVE_SUPPORT )
            list ( APPEND ADDITIONAL_NVCC_RELEASE_FLAGS -Xcompiler -march=native )
        endif ( MARCH_NATIVE_SUPPORT )
        
    endif ( WIN32 )
    
    set ( ADDITIONAL_NVCC_FLAGS "${ADDITIONAL_NVCC_FLAGS}" CACHE STRING "additional nvcc compiler flags" )
    set ( ADDITIONAL_NVCC_RELEASE_FLAGS "${ADDITIONAL_NVCC_RELEASE_FLAGS}" CACHE STRING "additional nvcc release compiler flags" )
    mark_as_advanced ( ADDITIONAL_NVCC_FLAGS ADDITIONAL_NVCC_RELEASE_FLAGS )
        
    list ( APPEND CUDA_NVCC_FLAGS "${ADDITIONAL_NVCC_FLAGS}" )
    list ( APPEND CUDA_NVCC_FLAGS_RELEASE "${ADDITIONAL_NVCC_RELEASE_FLAGS}" )
    
    # TODO: determine cuda compute capability and use highest
    # with sm_20 no warnings about Cannot tell what pointer points to, assuming global memory space in Release build
    # We need at least compute capability 1.3, so if no architecture is specified set it here
    if ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    	list ( APPEND CUDA_NVCC_FLAGS "-arch=sm_${CUDA_COMPUTE_CAPABILITY}" )
    endif ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    
endif( CUDA_FOUND AND LAMA_USE_CUDA )
