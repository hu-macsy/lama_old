#
# Find the hmemo (hybrid memory) includes and libraries
#
# SCAI_HMEMO_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_HMEMO_INCLUDE_DIR - the memory include dir
# SCAI_HMEMO_LIBRARY     - libraries to link against

if ( NOT SCAI_HMEMO_INCLUDE_DIR )
    find_path ( SCAI_HMEMO_INCLUDE_DIR hmemo.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_HMEMO_INCLUDE_PATH}/scai
        ${SCAI_HMEMO_ROOT}/include/scai
    )
    get_filename_component ( SCAI_HMEMO_INCLUDE_DIR ${SCAI_HMEMO_INCLUDE_DIR} PATH )
endif ( NOT SCAI_HMEMO_INCLUDE_DIR )

set ( SCAI_HMEMO_INCLUDE_DIR ${SCAI_HMEMO_INCLUDE_DIR} CACHE PATH "Path to HMEMO include dir" FORCE )

find_library ( SCAI_HMEMO_LIBRARY scai_hmemo
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_HMEMO_LIBRARY_PATH}
    ${SCAI_HMEMO_ROOT}/lib
)

set ( SCAI_HMEMO_FOUND FALSE )
if ( SCAI_HMEMO_INCLUDE_DIR )
    if ( SCAI_HMEMO_LIBRARY)
        set ( SCAI_HMEMO_FOUND TRUE )
    endif ( SCAI_HMEMO_LIBRARY )
endif ( SCAI_HMEMO_INCLUDE_DIR)

mark_as_advanced ( SCAI_HMEMO_FOUND SCAI_HMEMO_INCLUDE_DIR SCAI_HMEMO_LIBRARY )
