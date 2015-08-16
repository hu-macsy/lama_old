#
# Find the memory includes and libraries
#
# SCAI_MEMORY_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_MEMORY_INCLUDE_DIR - the memory include dir
# SCAI_MEMORY_LIBRARY     - libraries to link against

if    ( NOT DEFINED SCAI_MEMORY_INCLUDE_DIR )
    find_path ( SCAI_MEMORY_INCLUDE_DIR memory.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_MEMORY_INCLUDE_PATH}/scai
        ${SCAI_MEMORY_ROOT}/include/scai
    )
endif ( NOT DEFINED SCAI_MEMORY_INCLUDE_DIR )

find_library ( SCAI_MEMORY_LIBRARY scai_memory
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_MEMORY_LIBRARY_PATH}
    ${SCAI_MEMORY_ROOT}/lib
)

if    ( SCAI_MEMORY_INCLUDE_DIR )
    if    (SCAI_MEMORY_LIBRARY)
        set ( SCAI_MEMORY_FOUND TRUE )
    endif ( SCAI_MEMORY_LIBRARY )
endif (SCAI_MEMORY_INCLUDE_DIR)

mark_as_advanced ( SCAI_MEMORY_FOUND SCAI_MEMORY_INCLUDE_DIR SCAI_MEMORY_LIBRARY )
