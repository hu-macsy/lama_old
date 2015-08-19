#
# Find the logging includes and libraries
#
# SCAI_LOGGING_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_LOGGING_INCLUDE_DIR - the logging include dir
# SCAI_LOGGING_LIBRARY     - libraries to link against
# SCAI_LOGGING_LEVEL       - level of logging, e.g.TRACE, DEBUG, INFO
# SCAI_LOGGING_FLAG        - compile flag for logging

if ( NOT SCAI_LOGGING_INCLUDE_DIR )
    find_path ( SCAI_LOGGING_INCLUDE_DIR logging.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LOGGING_INCLUDE_PATH}/scai
        ${SCAI_LOGGING_ROOT}/include/scai
    )
endif ( NOT SCAI_LOGGING_INCLUDE_DIR )

set ( SCAI_LOGGING_INCLUDE_DIR ${SCAI_LOGGING_INCLUDE_DIR} CACHE PATH "Path to LOGGING include dir" FORCE )

find_library ( SCAI_LOGGING_LIBRARY scai_logging
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_LOGGING_LIBRARY_PATH}
    ${SCAI_LOGGING_ROOT}/lib
)

if ( SCAI_LOGGING_INCLUDE_DIR )
    if ( SCAI_LOGGING_LIBRARY)
        set ( SCAI_LOGGING_FOUND TRUE )
    endif ( SCAI_LOGGING_LIBRARY )
endif ( SCAI_LOGGING_INCLUDE_DIR)

# Using the logging library requires setting of SCAI_LOG_LEVEL_xxx at compile time
# CMake variable SCAI_LOGGING_LEVEL (cache) is used for the correct setting

set( LOG_CHOICES "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "OFF" )

## LOGGING Level
#
#  Debug   : use -DSCAI_LOG_LEVEL_DEBUG
#  Release : use -DSCAI_LOG_LEVEL_INFO
#  
#  For serious problems: -DLOG_LEVEL_TRACE
#  For benchmarks:       -DLOG_LEVEL_OFF (or -DLOG_LEVEL_FATAL, -DLOG_LEVEL_ERROR)

if ( NOT SCAI_LOGGING_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_LOG_LEVEL "INFO" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    else ( )
        set ( DEFAULT_LOG_LEVEL "TRACE" )
    endif ( )

    if ( NOT scai_logging_FIND_QUIETLY )
        message ( STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE} implies default logging level ${DEFAULT_LOG_LEVEL}" )
    endif ()

endif ( )

set ( SCAI_LOGGING_LEVEL ${DEFAULT_LOG_LEVEL} CACHE STRING
      "Choose level of compile time logging: ${LOG_CHOICES}" )

# check for a correct entry

set ( _okay FALSE )
foreach ( ITEM ${LOG_CHOICES} )
   if ( ${SCAI_LOGGING_LEVEL} MATCHES ${ITEM} )
      set ( _okay TRUE )
   endif ()
endforeach()

if ( NOT _okay )
    message ( FATAL_ERROR "Logging value ${SCAI_LOGGING_LEVEL} is invalid, use one of ${LOG_CHOICES}" )
endif ( NOT _okay )

set ( CACHE LOG_LEVEL PROPERTY STRINGS ${LOG_CHOICES} )

# compile flag for logging should not be put in the cache, avoids errors for wrong settings

set ( SCAI_LOGGING_FLAG "SCAI_LOG_LEVEL_${SCAI_LOGGING_LEVEL}" )

# Note: files using logging should be compiles with the corresponding flag
# add_definitions ( -D${SCAI_LOGGING_FLAG} )

mark_as_advanced ( SCAI_LOGGING_INCLUDE_DIR SCAI_LOGGING_LIBRARY )

