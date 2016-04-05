#
# Find the logging includes and libraries
#
# InputVariables:
# 
# CMAKE_INSTALL_PREFIX : directory is used to find the logging installation
# LOGGING_ROOT         : installation directory where the logging library is installed 
#
# SCAI_LOGGING_INCLUDE_PATH : environment variable used to find logging include directory
# SCAI_LOGGING_LIBRARY_PATH : environment variable used to find logging include directory
#
# OutputVariables:
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
    get_filename_component ( SCAI_LOGGING_INCLUDE_DIR ${SCAI_LOGGING_INCLUDE_DIR} PATH )
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

mark_as_advanced ( SCAI_LOGGING_INCLUDE_DIR SCAI_LOGGING_LIBRARY )

include ( Settings/logging )
