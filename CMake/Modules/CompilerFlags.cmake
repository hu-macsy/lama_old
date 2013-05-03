###
 # @file CompilerFlags.cmake
 #
 # @license
 # Copyright (c) 2013
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
###

#### needed to find BLAS Libraries ####

if ( NOT WIN32 )
    enable_language ( Fortran )
endif ( NOT WIN32)

if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    find_program ( XIAR xiar )
    if ( XIAR )
        set ( CMAKE_AR "${XIAR}" )
    endif ( XIAR )
    mark_as_advanced ( XIAR )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

# Define Compile Flags for LAMA
#

message ( STATUS "${CMAKE_CXX_COMPILER_ID} compiler" )

include ( CheckCCompilerFlag )

#### compiler independent flag definition ####

# standard compiler flags
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} " )

# default warning flags
set ( CXX_WARNING_FLAGS "${CXX_WARNING_FLAGS} " )

# default release flags
set ( CXX_RELEASE_FLAGS "${CXX_RELEASE_FLAGS} " )

# default optimization flags (release build)
set ( ADDITIONAL_CXX_RELEASE_FLAGS "${ADDITIONAL_CXX_RELEASE_FLAGS} -O3 " )
if ( MARCH_NATIVE_SUPPORT )
    set ( ADDITIONAL_CXX_RELEASE_FLAGS "${ADDITIONAL_CXX_RELEASE_FLAGS} -march=native " )
endif ( MARCH_NATIVE_SUPPORT )


#### compiler dependent flag definition ####

# GNU
if ( CMAKE_COMPILER_IS_GNUCXX )
    set ( ADDITIONAL_CXX_FLAGS "-Wl,--no-as-needed " )
    set ( ADDITIONAL_CXX_WARNING_FLAGS "-Wextra -Wall " ) # -pedantic -std=c++98 " ) # -march=core02
    set ( ADDITIONAL_CXX_RELEASE_FLAGS "-ffast-math -msse4a " )
endif ( CMAKE_COMPILER_IS_GNUCXX )


# INTEL
if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    set ( ADDITIONAL_CXX_FLAGS "-fPIC -std=c++0x -shared-intel -wd1478 " ) # suppress warning 1478: deprecated auto_ptr
    set ( ADDITIONAL_CXX_WARNING_FLAGS "-w2 -Wall -Wcheck -Werror-all " ) # -Werror-all Warnings/Errors. No Remarks.
    set ( ADDITIONAL_CXX_RELEASE_FLAGS "-ipo -no-prec-div -xHost " )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )


# PGI
if ( CMAKE_CXX_COMPILER_ID MATCHES PGI )
    set ( ADDITIONAL_CXX_FLAGS "-fPIC -Kieee -Mipa=libc -DBOOST_HAS_THREADS " ) # -std=c++0x 
    # Disable warning 1097 to avoid warnings from openmpi headers with gcc specific attributes
    set ( ADDITIONAL_CXX_WARNING_FLAGS "--display_error_number --diag_suppress1097 " )
    set ( ADDITIONAL_CXX_RELEASE_FLAGS "-fast " )
endif ( CMAKE_CXX_COMPILER_ID MATCHES PGI )


#### concluding all defined compiler flags to CMAKE_..._FLAGS ####
set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS} ${ADDITIONAL_CXX_WARNING_FLAGS} ${COVERAGE_FLAGS}")

if ( WIN32 )
    #Disable warnings about insecure cstdlib functions of which secure alternatives are only available on windows
    set ( CXX_WARNING_FLAGS "-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS " )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXX_WARNING_FLAGS}" )
endif ( WIN32 )

set ( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} ${CXX_WARNING_FLAGS} " )
set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_RELEASE_FLAGS} ${ADDITIONAL_CXX_RELEASE_FLAGS} " )

#### CUDA specific compiler flags ####

if ( CUDA_FOUND AND LAMA_USE_CUDA )
    
    ### choosing the right compute capability
    ### we just start from version 1.3 ( 1.0 - 1.2 is not supported )
    LIST ( APPEND CC_CHOICES "13" "20" "21" "30" "35" )
    set ( CUDA_COMPUTE_CAPABILITY "13" CACHE STRING "CUDA compute capability (supported up from 13)" )
	set ( CACHE CUDA_COMPUTE_CAPABILITY PROPERTY STRINGS ${CC_CHOICES} )
    checkValue( ${CUDA_COMPUTE_CAPABILITY} "${CC_CHOICES}" )
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
        
        # Intel compiler
        if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
            set ( ADDITIONAL_NVCC_FLAGS "---compiler-bindir ${CMAKE_CXX_COMPILER};${ADDITIONAL_NVCC_FLAGS}" )  
        endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
        
        #-Xcompiler;-fno-inline is used because of compability issues of CUDA with gcc-4.4
        if ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
      	    set ( ADDITIONAL_NVCC_FLAGS "-g;-G;-Xcompiler;-fPIC;${ADDITIONAL_NVCC_FLAGS}" )
        else ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
       	    set ( ADDITIONAL_NVCC_FLAGS "-Xcompiler;-fPIC;${ADDITIONAL_NVCC_FLAGS}" )
        endif ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
        
        # set -march=core02,-mmmx,-msse,-msse2,-msse3,-mssse3,-msse4a flaggs here
        if ( MARCH_NATIVE_SUPPORT )
            set ( ADDITIONAL_NVCC_RELEASE_FLAGS "-O3;-use_fast_math;-Xcompiler;-ffast-math;-Xcompiler;-fno-inline;-Xcompiler;-march=native;${ADDITIONAL_NVCC_FLAGS}" )
        else ( MARCH_NATIVE_SUPPORT )
            set ( ADDITIONAL_NVCC_RELEASE_FLAGS "-O3;-use_fast_math;-Xcompiler;-ffast-math;-Xcompiler;-fno-inline;${ADDITIONAL_NVCC_FLAGS}" )
        endif ( MARCH_NATIVE_SUPPORT )
        
    endif ( WIN32 )
    
    # TODO: determine cuda compute capability and use highest
    # with sm_20 no warnings about Cannot tell what pointer points to, assuming global memory space in Release build
    
    # We need at least compute capability 1.3, so if no architecture is specified set it here
    if ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    	set ( ADDITIONAL_NVCC_FLAGS "-arch=sm_${CUDA_COMPUTE_CAPABILITY};${ADDITIONAL_NVCC_FLAGS}" )
    endif ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    
    list ( APPEND CUDA_NVCC_FLAGS "${ADDITIONAL_NVCC_FLAGS}" )
    list ( APPEND CUDA_NVCC_FLAGS_RELEASE "${ADDITIONAL_NVCC_RELEASE_FLAGS}" )
    
    mark_as_advanced ( ADDITIONAL_NVCC_FLAGS )
    mark_as_advanced ( ADDITIONAL_NVCC_RELEASE_FLAGS )
    
    message ( STATUS "CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}" )
    message ( STATUS "CUDA_NVCC_FLAGS_RELEASE ${CUDA_NVCC_FLAGS_RELEASE}" )
    
endif( CUDA_FOUND AND LAMA_USE_CUDA )

## add variables to cache

set ( ADDITIONAL_CXX_FLAGS "${ADDITIONAL_CXX_FLAGS}" CACHE STRING "additional cxx compiler flags" )
set ( ADDITIONAL_CXX_WARNING_FLAGS "${ADDITIONAL_CXX_WARNING_FLAGS}" CACHE STRING "additional cxx compiler flags concerning warnings" )
set ( ADDITIONAL_CXX_RELEASE_FLAGS "${ADDITIONAL_CXX_RELEASE_FLAGS}" CACHE STRING "addtional cxx compiler flags for release optimizations" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS ADDITIONAL_CXX_WARNING_FLAGS ADDITIONAL_CXX_RELEASE_FLAGS )

## hide flags (we do not use) from the default CMake screen 

set ( CMAKE_CXX_FLAGS_MINSIZEREL "${CMAKE_CXX_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
set ( CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}" CACHE INTERNAL "" )
set ( CMAKE_C_FLAGS_MINSIZEREL "${CMAKE_C_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
set ( CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO}" CACHE INTERNAL "" )
set ( CMAKE_Fortran_FLAGS_MINSIZEREL "${CMAKE_Fortran_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
set ( CMAKE_Fortran_FLAGS_RELWITHDEB "${CMAKE_Fortran_FLAGS_RELWITHDEB}" CACHE INTERNAL "" )
set ( CMAKE_EXE_LINKER_FLAGS_MINSIZEREL "${CMAKE_EXE_LINKER_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
set ( CMAKE_EXE_LINKER_FLAGS_RELWITHDEB "${CMAKE_EXE_LINKER_FLAGS_RELWITHDEB}" CACHE INTERNAL "" )
set ( CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL "${CMAKE_MODULE_LINKER_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
set ( CMAKE_MODULE_LINKER_FLAGS_RELWITHDEB "${CMAKE_MODULE_LINKER_FLAGS_RELWITHDEB}" CACHE INTERNAL "" )
set ( CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL "${CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
set ( CMAKE_SHARED_LINKER_FLAGS_RELWITHDEB "${CMAKE_SHARED_LINKER_FLAGS_RELWITHDEB}" CACHE INTERNAL "" )
set ( CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL "${CMAKE_SHARED_LINKER_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )

if ( CUDA_FOUND  )
    set ( CUDA_NVCC_FLAGS_MINSIZEREL "${CUDA_NVCC_FLAGS_MINSIZEREL}" CACHE INTERNAL "" )
    set ( CUDA_NVCC_FLAGS_RELWITHDEBINFO "${CUDA_NVCC_FLAGS_RELWITHDEBINFO}" CACHE INTERNAL "" )
    set ( CUDA_GENERATED_OUTPUT_DIR "${CUDA_GENERATED_OUTPUT_DIR}" CACHE INTERNAL "" )
    set ( CUDA_SDK_ROOT_DIR "$CUDA_SDK_ROOT_DIR" CACHE INTERNAL "" )
endif ( CUDA_FOUND  )

if ( MPI_FOUND )
    set ( MPIEXEC_POSTFLAGS "${MPIEXEC_POSTFLAGS}" CACHE INTERNAL "" )
    set ( MPIEXEC_PREFLAGS "${MPIEXEC_PREFLAGS}" CACHE INTERNAL "" )
endif ( MPI_FOUND)