# Check if verbose mode for CMAKE is selected
if ( DEFINED LAMA_CMAKE_VERBOSE AND LAMA_CMAKE_VERBOSE )
    #set ( LAMA_CMAKE_VERBOSE TRUE )
else ()
    #set ( LAMA_CMAKE_VERBOSE FALSE )
endif( DEFINED LAMA_CMAKE_VERBOSE AND LAMA_CMAKE_VERBOSE )

set ( CMAKE_SYSTEM_LIBRARY_PATH )
set ( LAMA_ROOT_DIR "${CMAKE_SOURCE_DIR}/.." )

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

get_property ( FIND_LIB64 GLOBAL PROPERTY FIND_LIBRARY_USE_LIB64_PATHS )
message ( STATUS "FindLib64: " ${FIND_LIB64} )


# Choose Default CMAKE_BUILD_TYPE
if ( NOT CMAKE_BUILD_TYPE )
  # Can be: (RelWithDebInfo)
  set ( CMAKE_BUILD_TYPE Debug CACHE STRING 
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif ( NOT CMAKE_BUILD_TYPE )



# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )


# Find required packages
set ( REQUIRED_PACKAGES_TO_FIND
        LAMA_BLAS
        #add required packages here
    )
    
# Find optional packages
set ( OPTIONAL_PACKAGES_TO_FIND
        OpenMP
        Doxygen
        Threads
        #add optional packages here
    )
    
#CUDA Only works with GCC on Linux
#TODO: This needs to be checked on windows
if ( CMAKE_COMPILER_IS_GNUCC )
    set ( OPTIONAL_PACKAGES_TO_FIND
          ${OPTIONAL_PACKAGES_TO_FIND}
          CUDA
          #CUDA needs to be found before OpenCL because CUDA_INCLUDE_DIRS is used as a hint to find the OpenCL Headers
          #OpenCL
    )
endif ( CMAKE_COMPILER_IS_GNUCC )




## BUILDTYPE

# Choose Default CMAKE_BUILD_TYPE
if ( NOT CMAKE_BUILD_TYPE )
  # Can be: (RelWithDebInfo)
  set ( CMAKE_BUILD_TYPE Debug CACHE STRING 
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

message ( STATUS "Build type is set to " ${CMAKE_BUILD_TYPE} )

## LOGGING Level
#
#  Debug   : use -DLOG_LEVEL_DEBUG
#  Release : use -DLOG_LEVEL_INFO
#  
#  For serious problems: -DLOG_LEVEL_TRACE
#  For benchmarks:       -DLOG_LEVEL_OFF (or -DLOG_LEVEL_FATAL, -DLOG_LEVEL_ERROR)

if ( NOT LAMA_LOG_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_LOG_LEVEL "INFO" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    else ()
        set ( DEFAULT_LOG_LEVEL "TRACE" )
    endif ()
endif ( NOT LAMA_LOG_LEVEL )

set ( LAMA_LOG_LEVEL ${DEFAULT_LOG_LEVEL} CACHE STRING
      "Choose level of compile time logging: TRACE, DEBUG, INFO, WARN, ERROR, OFF" )

add_definitions ( -DLAMA_LOG_LEVEL_${LAMA_LOG_LEVEL} )



## ASSERT Level
#
#  Debug   : use -DASSERT_LEVEL_DEBUG
#  Release : use -DASSERT_LEVEL_ERROR
#  
#  For benchmarks:       -DASSERT_LEVEL_OFF

if ( NOT LAMA_ASSERT_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_ASSERT_LEVEL "ERROR" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    else ()
        set ( DEFAULT_ASSERT_LEVEL "DEBUG" )
    endif ()
endif ( NOT LAMA_ASSERT_LEVEL )

set ( LAMA_ASSERT_LEVEL ${DEFAULT_ASSERT_LEVEL} CACHE STRING
      "Choose level of ASSERT: DEBUG, ERROR, OFF" )

add_definitions ( -DLAMA_ASSERT_LEVEL_${LAMA_ASSERT_LEVEL} )



## LAMA TRACE LEVEL
#
# If TRACE is set to OFF all LAMA_REGION macros in the code are
# completely ignored. If TRACE is set to VT, regions will be traced
# (entry, exit event) for VampirTrace.

if ( NOT LAMA_TRACE_LEVEL )
    set ( DEFAULT_TRACE_LEVEL "OFF" )
endif ( NOT LAMA_TRACE_LEVEL )

set ( LAMA_TRACE_LEVEL ${DEFAULT_TRACE_LEVEL} CACHE STRING
     "Choose level of TRACE: VT (for VampirTrace), TIME(region timing), SIMPLE(simple timing) or OFF (default)" )

add_definitions( -DLAMA_TRACE_LEVEL_${LAMA_TRACE_LEVEL} )


    