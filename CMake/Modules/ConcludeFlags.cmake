# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

# for static/dynamic linking
if    ( ${SCAI_LIBRARY_TYPE} MATCHES "STATIC" )
	set ( SCAI_START_LINK_LIBRARIES "-Wl,--whole-archive" )
	set ( SCAI_END_LINK_LIBRARIES "-Wl,--no-whole-archive" )
else  ( ${SCAI_LIBRARY_TYPE} MATCHES "STATIC" )
	set ( SCAI_START_LINK_LIBRARIES "-Wl,--no-as-needed" )
	set ( SCAI_END_LINK_LIBRARIES "-Wl,--as-needed" )
endif ( ${SCAI_LIBRARY_TYPE} MATCHES "STATIC" )

## add variables to cache with new names so they can be modified by the user via CCMAKE

# moved to packages
#set ( ADDITIONAL_CXX_FLAGS_LANG          "${SCAI_LANG_FLAGS}"          CACHE STRING "Addition language flags for using C++11 (if compiler capable)" )
#set ( ADDITIONAL_CXX_FLAGS_OPENMP        "${OpenMP_CXX_FLAGS}"         CACHE STRING "OpenMP flag (only if enabled)" )

set ( ADDITIONAL_CXX_FLAGS_CODE_COVERAGE "${SCAI_CODE_COVERAGE_FLAGS}" CACHE STRING "CXX flags used for code coverage (only if CC enabled)" )
set ( ADDITIONAL_CXX_FLAGS_DEBUG         "${SCAI_CXX_FLAGS_DEBUG}"     CACHE STRING "Addtional CXX compiler flags for Debug version" )
set ( ADDITIONAL_CXX_FLAGS_NO_OFFLOAD    "${MIC_NO_OFFLOAD_FLAG}"      CACHE STRING "MIC no offload flag (only if MIC disabled)" )
set ( ADDITIONAL_CXX_FLAGS_RELEASE       "${SCAI_CXX_FLAGS_RELEASE}"   CACHE STRING "Addtional CXX compiler flags for Release version" )
set ( ADDITIONAL_LINKER_FLAGS            "${SCAI_LINKER_FLAGS}"        CACHE STRING "Additional linker flags" )
set ( ADDITIONAL_WARNING_FLAGS           "${SCAI_WARNING_FLAGS}"       CACHE STRING "Compilation flags concerning warnings" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS_CODE_COVERAGE  ADDITIONAL_CXX_FLAGS_DEBUG  ADDITIONAL_CXX_FLAGS_RELEASE
                   ADDITIONAL_CXX_FLAGS_NO_OFFLOAD     #ADDITIONAL_CXX_FLAGS_LANG   ADDITIONAL_CXX_FLAGS_OPENMP
                   ADDITIONAL_LINKER_FLAGS             ADDITIONAL_WARNING_FLAGS
                 )

set ( CONCLUDE_CXX_FLAGS "${SCAI_CXX_FLAGS}" )

if    ( SCAI_COMMON_FOUND )
    set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${SCAI_COMMON_FLAGS}")
else  ( SCAI_COMMON_FOUND )

    if    ( USE_OPENMP )
        set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS_OPENMP}" )
    else  ( USE_OPENMP )
        set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS_NO_OPENMP}" )
    endif ( USE_OPENMP)

    if    ( CXX_SUPPORTS_C11 )
        set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS_LANG}" )
    endif ( CXX_SUPPORTS_C11 )

endif ( SCAI_COMMON_FOUND )

if    ( DEFINED USE_MIC )
    message( STATUS "DEFINED USE_MIC" )
    if    ( NOT USE_MIC )
        message( STATUS "NOT USE_MIC" )
        set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS_NO_OFFLOAD}" )
    endif ( NOT USE_MIC )
else  ( DEFINED USE_MIC )
    message( STATUS "NOT DEFINED USE_MIC" )
    set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS_NO_OFFLOAD}" )
endif ( DEFINED USE_MIC )

if    ( USE_CODE_COVERAGE )
    set ( CONCLUDE_CXX_FLAGS "${CONCLUDE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS_CODE_COVERAGE}")
endif ( USE_CODE_COVERAGE )

# remove leading and trailing whitespaces
string ( STRIP "${CONCLUDE_CXX_FLAGS}" CONCLUDE_CXX_FLAGS )

set ( CMAKE_CXX_FLAGS           "${CMAKE_CXX_FLAGS} ${CONCLUDE_CXX_FLAGS}" )
set ( CMAKE_CXX_FLAGS_RELEASE   "${CMAKE_CXX_FLAGS_RELEASE} ${ADDITIONAL_CXX_FLAGS_RELEASE} " )
set ( CMAKE_CXX_FLAGS_DEBUG     "${CMAKE_CXX_FLAGS_DEBUG} ${ADDITIONAL_CXX_FLAGS_DEBUG} " )
set ( CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )
set ( CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ADDITIONAL_LINKER_FLAGS} " )

# remove leading and trailing whitespaces
string ( STRIP "${CMAKE_CXX_FLAGS}"           CMAKE_CXX_FLAGS )
string ( STRIP "${CMAKE_CXX_FLAGS_DEBUG}"     CMAKE_CXX_FLAGS_DEBUG )
string ( STRIP "${CMAKE_CXX_FLAGS_RELEASE}"   CMAKE_CXX_FLAGS_RELEASE )
string ( STRIP "${CMAKE_EXE_LINKER_FLAGS}"    CMAKE_EXE_LINKER_FLAGS )
string ( STRIP "${CMAKE_SHARED_LINKER_FLAGS}" CMAKE_SHARED_LINKER_FLAGS )

message ( STATUS "### final CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}###" )

if ( CUDA_FOUND AND USE_CUDA )
    
    # TODO: determine cuda compute capability and use highest
    # with sm_20 no warnings about Cannot tell what pointer points to, assuming global memory space in Release build
    # We need at least compute capability 1.3, so if no architecture is specified set it here
    if    ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
        list ( APPEND CUDA_NVCC_FLAGS -arch=sm_${CUDA_COMPUTE_CAPABILITY} )
    endif ( NOT "${CUDA_NVCC_FLAGS}" MATCHES "-arch" )
    
    set ( ADDITIONAL_NVCC_FLAGS         "${SCAI_NVCC_FLAGS}"         CACHE STRING "additional nvcc compiler flags" )
    set ( ADDITIONAL_NVCC_FLAGS_DEBUG   "${SCAI_NVCC_FLAGS_DEBUG}"   CACHE STRING "additional nvcc debug compiler flags" )
    set ( ADDITIONAL_NVCC_FLAGS_RELEASE "${SCAI_NVCC_FLAGS_RELEASE}" CACHE STRING "additional nvcc release compiler flags" )
    mark_as_advanced ( ADDITIONAL_NVCC_FLAGS ADDITIONAL_NVCC_FLAGS_RELEASE ADDITIONAL_NVCC_FLAGS_DEBUG )

    list ( APPEND CUDA_NVCC_FLAGS         ${ADDITIONAL_NVCC_FLAGS} )
    list ( APPEND CUDA_NVCC_FLAGS_DEBUG   ${ADDITIONAL_NVCC_FLAGS_DEBUG} )
    list ( APPEND CUDA_NVCC_FLAGS_RELEASE ${ADDITIONAL_NVCC_FLAGS_RELEASE} )

    # remove leading and trailing whitespaces
    string ( STRIP "${CUDA_NVCC_FLAGS}"         CUDA_NVCC_FLAGS )
    string ( STRIP "${CUDA_NVCC_FLAGS_DEBUG}"   CUDA_NVCC_FLAGS_DEBUG )
    string ( STRIP "${CUDA_NVCC_FLAGS_RELEASE}" CUDA_NVCC_FLAGS_RELEASE )
    
endif ( CUDA_FOUND AND USE_CUDA )
