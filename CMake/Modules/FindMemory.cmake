#
# Find the memory includes and libraries
#
# MEMORY_FOUND       - Do not attempt to use if "no" or undefined
# MEMORY_INCLUDE_DIR - the memory include dir
# MEMORY_LIBRARY     - libraries to link against

if    ( NOT DEFINED MEMORY_INCLUDE_DIR )
    find_path ( MEMORY_INCLUDE_DIR memory.hpp
        /usr/local/include
        /usr/include
        ${CMAKE_INSTALL_PREFIX}/include
        $ENV{MEMORY_INCLUDE_PATH}
        ${MEMORY_ROOT}/include
    )
endif ( NOT DEFINED MEMORY_INCLUDE_DIR )

find_library ( MEMORY_LIBRARY memory
    /usr/local/lib
    /usr/lib
    $ENV{MEMORY_LIBRARY_PATH}
    ${MEMORY_ROOT}/lib
)

if    ( MEMORY_INCLUDE_DIR )
    if    (MEMORY_LIBRARY)
        set ( MEMORY_FOUND TRUE )
    endif ( MEMORY_LIBRARY )
endif (MEMORY_INCLUDE_DIR)

mark_as_advanced ( MEMORY_FOUND MEMORY_INCLUDE_DIR MEMORY_LIBRARY )