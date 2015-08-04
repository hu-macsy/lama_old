#
# Find the tasking includes and libraries
#
# TASKING_FOUND       - Do not attempt to use if "no" or undefined
# TASKING_INCLUDE_DIR - the tasking include dir
# TASKING_LIBRARY     - libraries to link against

if    ( NOT DEFINED TASKING_INCLUDE_DIR )
    find_path ( TASKING_INCLUDE_DIR tasking.hpp
        /usr/local/include
        /usr/include
        ${CMAKE_INSTALL_PREFIX}/include
        $ENV{TASKING_INCLUDE_PATH}
        ${TASKING_ROOT}/include
    )
endif ( NOT DEFINED TASKING_INCLUDE_DIR )

find_library ( TASKING_LIBRARY tasking
    /usr/local/lib
    /usr/lib
    $ENV{TASKING_LIBRARY_PATH}
    ${TASKING_ROOT}/lib
)

if    ( TASKING_INCLUDE_DIR )
    if    (TASKING_LIBRARY)
        set ( TASKING_FOUND TRUE )
    endif ( TASKING_LIBRARY )
endif (TASKING_INCLUDE_DIR)

mark_as_advanced ( TASKING_FOUND TASKING_INCLUDE_DIR TASKING_LIBRARY )