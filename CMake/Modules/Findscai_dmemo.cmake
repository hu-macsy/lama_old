#
# Find the dmemo (distributed memory) includes and libraries
#
# SCAI_DMEMO_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_DMEMO_INCLUDE_DIR - the memory include dir
# SCAI_DMEMO_LIBRARY     - libraries to link against

if ( NOT SCAI_DMEMO_INCLUDE_DIR )
    find_path ( SCAI_DMEMO_INCLUDE_DIR dmemo.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_DMEMO_INCLUDE_PATH}/scai
        ${SCAI_DMEMO_ROOT}/include/scai
    )
    get_filename_component ( SCAI_DMEMO_INCLUDE_DIR ${SCAI_DMEMO_INCLUDE_DIR} PATH )
endif ( NOT SCAI_DMEMO_INCLUDE_DIR )

set ( SCAI_DMEMO_INCLUDE_DIR ${SCAI_DMEMO_INCLUDE_DIR} CACHE PATH "Path to DMEMO include dir" FORCE )

find_library ( SCAI_DMEMO_LIBRARY scai_dmemo
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_DMEMO_LIBRARY_PATH}
    ${SCAI_DMEMO_ROOT}/lib
)

if ( SCAI_DMEMO_INCLUDE_DIR )
    if ( SCAI_DMEMO_LIBRARY)
        set ( SCAI_DMEMO_FOUND TRUE )
    endif ( SCAI_DMEMO_LIBRARY )
endif ( SCAI_DMEMO_INCLUDE_DIR)

mark_as_advanced ( SCAI_DMEMO_FOUND SCAI_DMEMO_INCLUDE_DIR SCAI_DMEMO_LIBRARY )
