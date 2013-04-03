
#### needed to find BLAS Libraries ###

if ( NOT WIN32 )
    enable_language ( Fortran )
endif ( NOT WIN32)

if ( CMAKE_C_COMPILER_ID MATCHES Intel )
    find_program ( XIAR xiar )
    if ( XIAR )
        set ( CMAKE_AR "${XIAR}" )
    endif ( XIAR )
    mark_as_advanced ( XIAR )
endif ( CMAKE_C_COMPILER_ID MATCHES Intel )





# Define Compile Flags for LAMA
#
# Variables which can be modified:
## [CXX|C]_WARNING_FLAGS         Warning flags for all CMAKE_BUILD_TYPEs
## [CXX|C]_RELEASE_FLAGS         Flags for release mode, compiler independent    
## ADDITIONAL_[C|CXX]_...        Compiler/Architecture specific flags

message ( STATUS "${CMAKE_CXX_COMPILER_ID} compiler" )

include ( CheckCCompilerFlag )

# default warning flags
set ( CXX_WARNING_FLAGS "" )
set ( C_WARNING_FLAGS "" )

# default optimization flags
#REMARK opt-flags are only for buidl type: release
set ( CXX_RELEASE_FLAGS " -O3 " )
set ( C_RELEASE_FLAGS " -O3 " )

if ( MARCH_NATIVE_SUPPORT )
    set ( CXX_RELEASE_FLAGS "-march=native " )
endif ( MARCH_NATIVE_SUPPORT )

# GNU

# gnu cxx
if ( CMAKE_COMPILER_IS_GNUCXX )
    set ( ADDITIONAL_CXX_WARNING_FLAGS "-Wextra -Wall -Werror -Wl,--no-as-needed" )#-pedantic -std=c++98 " )
    # -march=core02
    if ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL ppc64 )
        # compile in 32 Bit on Cell Broaband Engine Processor
        set ( ADDITIONAL_CXX_FLAGS "-m32 " )
        set ( CMAKE_REQUIRED_FLAGS "-m32 " )
        set ( ADDITIONAL_CXX_RELEASE_FLAGS "-ffast-math -maltivec " )
    else ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL ppc64 )
        set ( ADDITIONAL_CXX_RELEASE_FLAGS "-ffast-math -msse4a " )
    endif ( ${CMAKE_SYSTEM_PROCESSOR} STREQUAL ppc64 )
endif ( CMAKE_COMPILER_IS_GNUCXX )

# INTEL

# intel cxx
if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
    set ( ADDITIONAL_CXX_FLAGS "-std=c++0x -shared-intel " )
    set ( ADDITIONAL_CXX_WARNING_FLAGS "-w2 -Wall -Wcheck -Werror-all " ) # Warnings/Errors. No Remarks.
    set ( ADDITIONAL_CXX_RELEASE_FLAGS "-ipo -no-prec-div -xHost " )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

# PGI
if ( CMAKE_CXX_COMPILER_ID MATCHES PGI )
    # set ( ADDITIONAL_CXX_FLAGS "-std=c++0x " )
    # Disable warning 1097 to avoid warnings from openmpi headers with
    # gcc specific attributes
    set ( ADDITIONAL_CXX_WARNING_FLAGS "--display_error_number --diag_suppress1097 " )
    set ( ADDITIONAL_CXX_RELEASE_FLAGS "-fast " )
endif ( CMAKE_CXX_COMPILER_ID MATCHES PGI )


# profiling
if( CMAKE_PROFILE )
    if( CMAKE_COMPILER_IS_GNUCC )
        set ( ADDITIONAL_CXX_FLAGS "-pg " ${ADDITIONAL_CXX_FLAGS} )
    elseif ( CMAKE_C_COMPILER_ID MATCHES Intel )
        set ( ADDITIONAL_CXX_FLAGS "-p " ${ADDITIONAL_CXX_FLAGS} )
    endif ()
endif ( CMAKE_PROFILE )

# CONCLUDE
if ( NOT WIN32 )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS} ${ADDITIONAL_CXX_WARNING_FLAGS} ")
else ( NOT WIN32 )
    #Disable warnings about insecure cstdlib functions of which secure alternatives are only available on windows
    set ( CXX_WARNING_FLAGS "-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_WARNINGS " )
    set ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_CXX_FLAGS} ${CXX_WARNING_FLAGS}  ${ADDITIONAL_CXX_WARNING_FLAGS}" )
endif ( NOT WIN32 )
set ( CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${CXX_RELEASE_FLAGS} ${ADDITIONAL_CXX_RELEASE_FLAGS}" )



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

