#
# Find the common includes and libraries
#
# TRACING_FOUND       - Do not attempt to use if "no" or undefined
# TRACING_INCLUDE_DIR - the common include dir
# TRACING_LIBRARY     - libraries to link against

if    ( NOT DEFINED TRACING_INCLUDE_DIR )
    find_path ( TRACING_INCLUDE_DIR tracing.hpp
        /usr/local/include
        /usr/include
        ${CMAKE_INSTALL_PREFIX}/include
        $ENV{TRACING_INCLUDE_PATH}
        ${TRACING_ROOT}/include
    )
endif ( NOT DEFINED TRACING_INCLUDE_DIR )

find_library ( TRACING_LIBRARY tracing
    /usr/local/lib
    /usr/lib
    $ENV{TRACING_LIBRARY_PATH}
    ${TRACING_ROOT}/lib
)

if    ( TRACING_INCLUDE_DIR )
    if    (TRACING_LIBRARY)
        set ( TRACING_FOUND TRUE )
    endif ( TRACING_LIBRARY )
endif (TRACING_INCLUDE_DIR)

mark_as_advanced ( TRACING_FOUND TRACING_INCLUDE_DIR TRACING_LIBRARY )