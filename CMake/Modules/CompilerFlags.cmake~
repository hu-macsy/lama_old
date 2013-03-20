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
    set ( ADDITIONAL_CXX_WARNING_FLAGS "-Wextra -Wall -Wl,--no-as-needed" )#-pedantic -std=c++98 " )
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
