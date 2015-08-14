#
# Find the common includes and libraries
#
# SCAI_COMMON_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_COMMON_INCLUDE_DIR - the common include dir
# SCAI_COMMON_LIBRARY     - libraries to link against

if    ( NOT DEFINED SCAI_COMMON_INCLUDE_DIR )
    find_path ( SCAI_COMMON_INCLUDE_DIR common.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_COMMON_INCLUDE_PATH}/scai
        ${SCAI_COMMON_ROOT}/include/scai
    )
endif ( NOT DEFINED SCAI_COMMON_INCLUDE_DIR )

find_library ( SCAI_COMMON_LIBRARY scai_common
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_COMMON_LIBRARY_PATH}
    ${SCAI_COMMON_ROOT}/lib
)

if    ( SCAI_COMMON_INCLUDE_DIR )
    if    (SCAI_COMMON_LIBRARY)
        set ( SCAI_COMMON_FOUND TRUE )
    endif ( SCAI_COMMON_LIBRARY )
endif (SCAI_COMMON_INCLUDE_DIR)

mark_as_advanced ( SCAI_COMMON_FOUND SCAI_COMMON_INCLUDE_DIR SCAI_COMMON_LIBRARY )
