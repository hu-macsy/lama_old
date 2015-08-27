#### concluding all defined compiler flags to CMAKE_..._FLAGS ####

## scai common adds OpenMP_CXX_FLAGS and SCAI_LANG_FLAGS to SCAI_COMMON_FLAGS if found
if    ( SCAI_COMMON_FOUND ) 
	set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} ${SCAI_COMMON_FLAGS}" )
else  ( SCAI_COMMON_FOUND )
	set ( LAMA_CXX_FLAGS "${LAMA_CXX_FLAGS} ${OpenMP_CXX_FLAGS} ${SCAI_LANG_FLAGS}" )
endif ( SCAI_COMMON_FOUND ) 

## add variables to cache with new names so they can be modified by the user via CCMAKE

set ( ADDITIONAL_CXX_FLAGS "${LAMA_CXX_FLAGS}" CACHE STRING "additional flags for cxx compile and link" )
set ( ADDITIONAL_WARNING_FLAGS "${SCAI_WARNING_FLAGS}" CACHE STRING "compilation flags concerning warnings" )
set ( ADDITIONAL_CXX_FLAGS_RELEASE "${LAMA_CXX_FLAGS_RELEASE}" CACHE STRING "addtional cxx compiler flags for release optimizations" )
set ( ADDITIONAL_LINKER_FLAGS "${LAMA_LINKER_FLAGS}" CACHE STRING "additional linker flags" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS ADDITIONAL_WARNING_FLAGS ADDITIONAL_CXX_FLAGS_RELEASE ADDITIONAL_LINKER_FLAGS )

set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS}")
set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS_RELEASE} " )
set ( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )

# remove leading and trailing whitespaces
string ( STRIP "${CMAKE_CXX_FLAGS}" CMAKE_CXX_FLAGS )
string ( STRIP "${CMAKE_CXX_FLAGS_RELEASE}" CMAKE_CXX_FLAGS_RELEASE )
string ( STRIP "${CMAKE_EXE_LINKER_FLAGS}" CMAKE_EXE_LINKER_FLAGS )
string ( STRIP "${CMAKE_SHARED_LINKER_FLAGS}" CMAKE_SHARED_LINKER_FLAGS )

if ( CUDA_FOUND AND USE_CUDA )
    
    # TODO: determine cuda compute capability and use highest
    # with sm_20 no warnings about Cannot tell what pointer points to, assuming global memory space in Release build
    # We need at least compute capability 1.3, so if no architecture is specified set it here
    if    ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    	list ( APPEND CUDA_NVCC_FLAGS -arch=sm_${CUDA_COMPUTE_CAPABILITY} )
    endif ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    
    set ( ADDITIONAL_NVCC_FLAGS "${ADDITIONAL_NVCC_FLAGS}" CACHE STRING "additional nvcc compiler flags" )
    set ( ADDITIONAL_NVCC_RELEASE_FLAGS "${ADDITIONAL_NVCC_RELEASE_FLAGS}" CACHE STRING "additional nvcc release compiler flags" )
    mark_as_advanced ( ADDITIONAL_NVCC_FLAGS ADDITIONAL_NVCC_RELEASE_FLAGS )

    list ( APPEND CUDA_NVCC_FLAGS ${ADDITIONAL_NVCC_FLAGS} )
    list ( APPEND CUDA_NVCC_FLAGS_RELEASE ${ADDITIONAL_NVCC_RELEASE_FLAGS} )
    
    # remove leading and trailing whitespaces
    string ( STRIP "${CUDA_NVCC_FLAGS}" CUDA_NVCC_FLAGS )
    string ( STRIP "${CUDA_NVCC_FLAGS_RELEASE}" CUDA_NVCC_FLAGS_RELEASE )
    
endif ( CUDA_FOUND AND USE_CUDA )