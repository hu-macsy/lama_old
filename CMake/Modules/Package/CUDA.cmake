###
 # @file CUDA.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief All search and settings for CUDA concentrated
 # @author Lauretta Schubert
 # @date 03.08.2015
###

### CUDA_FOUND            - if CUDA is found
### USE_CUDA              - if CUDA is enabled
### CUDA_ENABLED          - if CUDA_FOUND and USE_CUDA
### SCAI_CUDA_INCLUDE_DIR - CUDA include directory
### SCAI_CUDA_LIBRARIES   - all needed CUDA libraries

set ( CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Host side compiler used by NVCC" )

find_package ( CUDA ${SCAI_FIND_PACKAGE_FLAGS} )

# find out if host compiler version is supported by CUDA installation
include ( Compiler/cuda/CheckHostCompilerCompatibility )

# find out CUDA compute capability by test program
include ( Compiler/cuda/ComputeCapabilityCheck )

# ALLOW to switch off CUDA explicitly
include ( Functions/setAndCheckCache )
setAndCheckCache ( CUDA )
set ( USE_CUDA ${USE_CUDA} CACHE BOOL "Enable / Disable use of CUDA" )

set ( CUDA_ENABLED FALSE )
if ( CUDA_FOUND AND USE_CUDA )
	set ( CUDA_ENABLED TRUE )

	if    ( ${CUDA_VERSION_MAJOR} LESS 5 )
		message ( FATAL_ERROR "LAMA supports CUDA with SDK greater 5.0, your installation is ${CUDA_VERSION}. Use a newer CUDA installation or disable CUDA." )
	endif ( ${CUDA_VERSION_MAJOR} LESS 5 )

	# set nvcc compiler flags, if not added as external project (has flags from parent)
	if    ( NOT SCAI_COMPLETE_BUILD )
		include ( Compiler/cuda/SetNVCCFlags )
	endif ( NOT SCAI_COMPLETE_BUILD )
	
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

# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( CUDA_TOOLKIT_ROOT_DIR CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD CUDA_HOST_COMPILER )

# LAMA irrelevant entries will be removed from cmake GUI completely
set ( CUDA_BUILD_CUBIN "${CUDA_BUILD_CUBIN}" CACHE INTERNAL "" )
set ( CUDA_BUILD_EMULATION "${CUDA_BUILD_EMULATION}" CACHE INTERNAL "" )
set ( CUDA_SEPARABLE_COMPILATION "${CUDA_SEPARABLE_COMPILATION}" CACHE INTERNAL "" )
# CUDA_64_BIT_DEVICE_CODE CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE CUDA_GENERATED_OUTPUT_DIR CUDA_TARGET_CPU_ARCH
# CUDA_TOOLKIT_INCLUDE CUDA_TOOLKIT_ROOT_DIR CUDA_TOOLKIT_TARGET_DIR
# CUDA_cublasemu_LIBRARY CUDA_cufft_LIBRARY CUDA_cufftemu_LIBRARY

if    ( USE_CUDA AND NOT CUDA_FOUND )
    message( FATAL_ERROR "Build of LAMA Cuda enabled, but configuration is incomplete!")
endif ( USE_CUDA AND NOT CUDA_FOUND )
