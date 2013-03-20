# Build a CUDA library that will be linked with lama++

### set include directories ### 

include_directories ( . ${CUDA_INCLUDE_DIRS} )


### compiler flags ###

if ( CUDA_FOUND )
    set ( CUDA_VERBOSE_BUILD OFF )
    set ( CUDA_BUILD_EMULATION OFF )
    # unfortunately we can not propagate the host flags to CUDA
    # because this issues to much warning in cuda headers
    # TODO: maybe we can change this with future CUDA releases
    if ( WIN32 )
		set ( CUDA_PROPAGATE_HOST_FLAGS ON )
		set ( CUDA_NVCC_FLAGS "-arch=sm_13" )
    else ( WIN32 )
		set ( CUDA_PROPAGATE_HOST_FLAGS OFF )
		#-Xcompiler;-fno-inline is used because of compability issues of CUDA with gcc-4.4
		if ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
		    set ( CUDA_NVCC_FLAGS "-arch=sm_13;-g;-G;-Xcompiler;-fPIC" )
		else ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
		    set ( CUDA_NVCC_FLAGS "-arch=sm_13;-Xcompiler;-fPIC" )
		endif ( ${CMAKE_BUILD_TYPE} MATCHES "Debug" )
		# set -march=core02,-mmmx,-msse,-msse2,-msse3,-mssse3,-msse4a flaggs here
		if ( MARCH_NATIVE_SUPPORT )
			set ( CUDA_NVCC_FLAGS_RELEASE "-O3;-use_fast_math;-Xcompiler;-ffast-math;-Xcompiler;-fno-inline;-Xcompiler;-march=native" )
		else ( MARCH_NATIVE_SUPPORT )
			set ( CUDA_NVCC_FLAGS_RELEASE "-O3;-use_fast_math;-Xcompiler;-ffast-math;-Xcompiler;-fno-inline" )
		endif ( MARCH_NATIVE_SUPPORT )
    endif ( WIN32 )

endif( CUDA_FOUND )

### set source files ###

set ( CXXSOURCES
        CUDAInterface.cpp
        CUDAStreamSyncToken.cpp
        CUDAContext.cpp
        CUDAContextManager.cpp
        CUDAError.cpp
        CUDAHostContext.cpp
        CUDAHostContextManager.cpp
        CUDAException.cpp
        CUDABLAS1.cpp
        CUDABLAS2.cpp
        CUDABLAS3.cpp
        CUDALAPACK.cpp
        CUDATracerSyncToken.cpp
    )
    
set ( CUDA_SOURCES
        blas1.cu
        blas3.cu
        emptyKernel.cu
        CUDAUtils.cu
        CUDACSRUtils.cu
        CUDAELLUtils.cu
        CUDAJDSUtils.cu
        CUDADIAUtils.cu
        CUDACOOUtils.cu
    )

### cuda compile ###

if ( NOT WIN32 )
    cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
endif( NOT WIN32 )

### add library ###

add_library( lamacuda SHARED ${CXXSOURCES} ${CUDA_FILES} )

### target link libraries ###

target_link_libraries( lamacuda ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse.so ${CUDA_cublas_LIBRARY} ${CUDA_CUDART_LIBRARY} ${CUDA_CUDA_LIBRARY} )

### install library

install ( TARGETS lamacuda DESTINATION lib )

