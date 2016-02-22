#
# Find the lama includes and libraries
#
# SCAI_LAMA_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_LAMA_INCLUDE_DIR - the memory include dir
# SCAI_LAMA_LIBRARY     - libraries to link against

if ( NOT SCAI_LAMA_INCLUDE_DIR )
    find_path ( SCAI_LAMA_INCLUDE_DIR lama.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LAMA_INCLUDE_PATH}/scai
        ${SCAI_LAMA_ROOT}/include/scai
    )
endif ( NOT SCAI_LAMA_INCLUDE_DIR )

set ( SCAI_LAMA_INCLUDE_DIR ${SCAI_LAMA_INCLUDE_DIR} CACHE PATH "Path to LAMA include dir" FORCE )

find_library ( SCAI_LAMA_LIBRARY scai_lama
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_LAMA_LIBRARY_PATH}
    ${SCAI_LAMA_ROOT}/lib
)

if ( SCAI_LAMA_INCLUDE_DIR )
    if ( SCAI_LAMA_LIBRARY)
        set ( SCAI_LAMA_FOUND TRUE )
    endif ( SCAI_LAMA_LIBRARY )
endif ( SCAI_LAMA_INCLUDE_DIR)

mark_as_advanced ( SCAI_LAMA_FOUND SCAI_LAMA_INCLUDE_DIR SCAI_LAMA_LIBRARY )
