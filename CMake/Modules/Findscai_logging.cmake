#
# Find the logging includes and libraries
#
# SCAI_LOGGING_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_LOGGING_INCLUDE_DIR - the logging include dir
# SCAI_LOGGING_LIBRARY     - libraries to link against

if    ( NOT DEFINED SCAI_LOGGING_INCLUDE_DIR )
    find_path ( SCAI_LOGGING_INCLUDE_DIR logging.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LOGGING_INCLUDE_PATH}/scai
        ${SCAI_LOGGING_ROOT}/include/scai
    )
endif ( NOT DEFINED SCAI_LOGGING_INCLUDE_DIR )

find_library ( SCAI_LOGGING_LIBRARY scai_logging
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_LOGGING_LIBRARY_PATH}
    ${SCAI_LOGGING_ROOT}/lib
)

if    ( SCAI_LOGGING_INCLUDE_DIR )
    if    (SCAI_LOGGING_LIBRARY)
        set ( SCAI_LOGGING_FOUND TRUE )
    endif ( SCAI_LOGGING_LIBRARY )
endif (SCAI_LOGGING_INCLUDE_DIR)

mark_as_advanced ( SCAI_LOGGING_FOUND SCAI_LOGGING_INCLUDE_DIR SCAI_LOGGING_LIBRARY )
