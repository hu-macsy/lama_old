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
 # @brief All search and settings for CUDA concentrated
 # @author Lauretta Schubert
 # @date 03.08.2015
 # @since 2.0.0
###

### CUDA_FOUND            - if CUDA is found
### USE_CUDA              - if CUDA is enabled
### SCAI_CUDA_INCLUDE_DIR - CUDA include directory
### SCAI_CUDA_LIBRARIES   - all needed CUDA libraries

set ( CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Host side compiler used by NVCC" )

find_package ( CUDA ${SCAI_FIND_PACKAGE_FLAGS} )

# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( CUDA_TOOLKIT_ROOT_DIR CUDA_BUILD_CUBIN CUDA_BUILD_EMULATION CUDA_SDK_ROOT_DIR
				   CUDA_VERBOSE_BUILD CUDA_HOST_COMPILER CUDA_SEPARABLE_COMPILATION )

# ALLOW to switch off CUDA explicitly
include ( Functions/setAndCheckCache )
setAndCheckCache ( CUDA )

if ( CUDA_FOUND AND USE_CUDA )

	# find out CUDA compute capability by test program
	include ( CUDAComputeCapability )
	
	# set nvcc compiler flags, if not added as external project (has flags from parent)
	if    ( NOT SCAI_COMPLETE_BUILD )
		include ( SetNVCCFlags )
	endif ( NOT SCAI_COMPLETE_BUILD )
	
	### Check for cuSPASE library, Version 2 (since CUDA 5.0)
	
	if ( CUDA_VERSION_MAJOR MATCHES "5" )
	
	    #message( STATUS "Check for cuSPARSE V2 include file in ${CUDA_INCLUDE_DIRS}" )
	    
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
	
	### Older cmake version have not set CUDA_cusparse_LIBRARY
	
	if ( NOT CUDA_cusparse_LIBRARY )
	
	    ### cusparse is usually in same directory as cublas
	
	    get_filename_component( HINT_CUDA_LIBRARY_DIR ${CUDA_cublas_LIBRARY} PATH )
	
	    find_library( CUDA_cusparse_LIBRARY NAMES cusparse
	                  HINTS ${HINT_CUDA_LIBRARY_DIR} )
	
	    mark_as_advanced( CUDA_cusparse_LIBRARY )
	
	endif ( NOT CUDA_cusparse_LIBRARY )
	
	# just for making it the same variable ending for all packages
	set ( SCAI_CUDA_INCLUDE_DIR ${CUDA_INCLUDE_DIRS} )
	
	# conclude all needed CUDA libraries
	set ( SCAI_CUDA_LIBRARIES ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} )

	get_filename_component ( SCAI_CUDA_LIBRARY_PATH ${CUDA_CUDART_LIBRARY} PATH CACHE )

endif ( CUDA_FOUND AND USE_CUDA )

if    ( USE_CUDA AND NOT CUDA_FOUND )
    message( FATAL_ERROR "Build of LAMA Cuda enabled, but configuration is incomplete!")
endif ( USE_CUDA AND NOT CUDA_FOUND )
