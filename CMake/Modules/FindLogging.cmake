#
# Find the logging includes and libraries
#
# LOGGING_FOUND       - Do not attempt to use if "no" or undefined
# LOGGING_INCLUDE_DIR - the logging include dir
# LOGGING_LIBRARY     - libraries to link against

if    ( NOT DEFINED LOGGING_INCLUDE_DIR )
    find_path ( LOGGING_INCLUDE_DIR logging.hpp
        /usr/local/include
        /usr/include
        ${CMAKE_INSTALL_PREFIX}/include
        $ENV{LOGGING_INCLUDE_PATH}
        ${LOGGING_ROOT}/include
    )
endif ( NOT DEFINED LOGGING_INCLUDE_DIR )

find_library ( LOGGING_LIBRARY logging
    /usr/local/lib
    /usr/lib
    $ENV{LOGGING_LIBRARY_PATH}
    ${LOGGING_ROOT}/lib
)

if    ( LOGGING_INCLUDE_DIR )
    if    (LOGGING_LIBRARY)
        set ( LOGGING_FOUND TRUE )
    endif ( LOGGING_LIBRARY )
endif (LOGGING_INCLUDE_DIR)

mark_as_advanced ( LOGGING_FOUND LOGGING_INCLUDE_DIR LOGGING_LIBRARY )