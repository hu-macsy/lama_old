#
# Find the common includes and libraries
#
# COMMON_FOUND       - Do not attempt to use if "no" or undefined
# COMMON_INCLUDE_DIR - the common include dir
# COMMON_LIBRARY     - libraries to link against

if    ( NOT DEFINED COMMON_INCLUDE_DIR )
    find_path ( COMMON_INCLUDE_DIR common.hpp
        /usr/local/include
        /usr/include
        ${CMAKE_INSTALL_PREFIX}/include
        $ENV{COMMON_INCLUDE_PATH}
        ${COMMON_ROOT}/include
    )
endif ( NOT DEFINED COMMON_INCLUDE_DIR )

find_library ( COMMON_LIBRARY common
    /usr/local/lib
    /usr/lib
    $ENV{COMMON_LIBRARY_PATH}
    ${COMMON_ROOT}/lib
)

if    ( COMMON_INCLUDE_DIR )
    if    (COMMON_LIBRARY)
        set ( COMMON_FOUND TRUE )
    endif ( COMMON_LIBRARY )
endif (COMMON_INCLUDE_DIR)

mark_as_advanced ( COMMON_FOUND COMMON_INCLUDE_DIR COMMON_LIBRARY )