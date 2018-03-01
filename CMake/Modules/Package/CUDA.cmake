###
 # @file CUDA.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief All search and settings for CUDA concentrated
 # @author Lauretta Schubert
 # @date 03.08.2015
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

### CUDA_FOUND            - if CUDA is found
### USE_CUDA              - if CUDA is enabled
### SCAI_CUDA_INCLUDE_DIR - CUDA include directory
### SCAI_CUDA_LIBRARIES   - all needed CUDA libraries

### As several SCAI modules might depend on CUDA we try to prevent
### callling this configuration file several times

scai_pragma_once ()

set ( CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE FILEPATH "Host side compiler used by NVCC" )

# MINIMUM_VERSION 7.0 because we require support of -std=c++11 

find_package ( CUDA ${SCAI_FIND_PACKAGE_FLAGS} 7.0 )

# find out if host compiler version is supported by CUDA installation
include ( Compiler/cuda/CheckHostCompilerCompatibility )

# find out CUDA compute capability by test program
include ( Compiler/cuda/ComputeCapabilityCheck )

# ALLOW to switch off CUDA explicitly

scai_build_variable ( NAME      USE_CUDA
                      BOOL 
                      DEFAULT   ${CUDA_FOUND}
                      DOCSTRING "use of CUDA (for NVIDIA GPUs)" )

# define own FIND_CUDA_HELPER_LIBS macro as it is use in cmakes FindCUDA module
macro( SCAI_FIND_CUDA_HELPER_LIBS _name) # rename so it is not confused
    cuda_find_library_local_first(CUDA_${_name}_LIBRARY ${_name} "\"${_name}\" library")
    mark_as_advanced(CUDA_${_name}_LIBRARY)
endmacro()

if ( CUDA_FOUND AND USE_CUDA )

    if ( ${CUDA_VERSION_MAJOR} LESS 7 )
        message ( FATAL_ERROR "LAMA supports CUDA with SDK at least 7.0, your installation is ${CUDA_VERSION}. Use a newer CUDA installation or disable CUDA." )
    endif ()

    # set nvcc compiler flags, if not added as external project (has flags from parent)
    if ( NOT SCAI_COMPLETE_BUILD )
        include ( Compiler/cuda/SetNVCCFlags )
    endif ()
    
    ### find cuda helper libraries when cmake does not
    # get directory ov blas library as hint
    get_filename_component( HINT_CUDA_LIBRARY_DIR ${CUDA_cublas_LIBRARY} PATH )

    # cusparse showed up in version 3.2, so it should be there
    if ( NOT CUDA_cusparse_LIBRARY )
        ### cusparse is usually in same directory as cublas
        find_library( CUDA_cusparse_LIBRARY NAMES cusparse HINTS ${HINT_CUDA_LIBRARY_DIR} )
        mark_as_advanced( CUDA_cusparse_LIBRARY )
    endif ( NOT CUDA_cusparse_LIBRARY )

    #endif ( CUDA_VERSION VERSION_LESS "3.2" )
    
    # cusolver showed up in version 7.0 that is now already minimum
    if ( NOT CUDA_cusolver_LIBRARY )
        ### cusolver is usually in same directory as cublas
        find_library( CUDA_cusolver_LIBRARY NAMES cusolver HINTS ${HINT_CUDA_LIBRARY_DIR} )
        mark_as_advanced( CUDA_cusolver_LIBRARY )
    endif ( NOT CUDA_cusolver_LIBRARY ) 

    # just for making it the same variable ending for all packages
    set ( SCAI_CUDA_INCLUDE_DIR ${CUDA_INCLUDE_DIRS} )
    
    # conclude all needed CUDA libraries
    set ( SCAI_CUDA_LIBRARIES ${CUDA_CUDA_LIBRARY} ${CUDA_CUDART_LIBRARY} 
                              ${CUDA_cublas_LIBRARY} ${CUDA_cusparse_LIBRARY} 
                              ${CUDA_cufft_LIBRARY} ${CUDA_cusolver_LIBRARY} )

    get_filename_component ( SCAI_CUDA_LIBRARY_PATH ${CUDA_CUDART_LIBRARY} PATH CACHE )

endif ( CUDA_FOUND AND USE_CUDA )

# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( CUDA_TOOLKIT_ROOT_DIR CUDA_SDK_ROOT_DIR CUDA_VERBOSE_BUILD CUDA_HOST_COMPILER )

# LAMA irrelevant entries will be removed from cmake GUI completely

set ( CUDA_BUILD_CUBIN "${CUDA_BUILD_CUBIN}" CACHE INTERNAL "" )
set ( CUDA_BUILD_EMULATION "${CUDA_BUILD_EMULATION}" CACHE INTERNAL "" )
set ( CUDA_SEPARABLE_COMPILATION "${CUDA_SEPARABLE_COMPILATION}" CACHE INTERNAL "" )
set ( CUDA_USE_STATIC_CUDA_RUNTIME "${CUDA_USE_STATIC_CUDA_RUNTIME}" CACHE INTERNAL "" )
set ( CUDA_rt_LIBRARY "${CUDA_rt_LIBRARY}" CACHE INTERNAL "" )

if ( USE_CUDA AND NOT CUDA_FOUND )
    message( FATAL_ERROR "Build of LAMA Cuda enabled, but configuration is incomplete!")
endif ()

scai_summary_external ( NAME       CUDA
                        ENABLED    ${USE_CUDA}
                        FOUND      ${CUDA_FOUND} 
                        VERSION    ${CUDA_VERSION} 
                        INCLUDE    ${SCAI_CUDA_INCLUDE_DIR} 
                        LIBRARIES  ${SCAI_CUDA_LIBRARIES} 
                        CXX_FLAGS  "optimize for compute capability ${CUDA_COMPUTE_CAPABILITY}"
                        EXECUTABLE "${CUDA_NVCC_EXECUTABLE}"
                      )
