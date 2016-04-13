#
# Find the utilskernel includes and libraries
#
# SCAI_UTILSKERNEL_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_UTILSKERNEL_INCLUDE_DIR - the memory include dir
# SCAI_UTILSKERNEL_LIBRARY     - libraries to link against

if ( NOT SCAI_UTILSKERNEL_INCLUDE_DIR )
    find_path ( SCAI_UTILSKERNEL_INCLUDE_DIR utilskernel.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_UTILSKERNEL_INCLUDE_PATH}/scai
        ${SCAI_UTILSKERNEL_ROOT}/include/scai
    )
endif ( NOT SCAI_UTILSKERNEL_INCLUDE_DIR )

set ( SCAI_UTILSKERNEL_INCLUDE_DIR ${SCAI_UTILSKERNEL_INCLUDE_DIR} CACHE PATH "Path to UTILSKERNEL include dir" FORCE )

find_library ( SCAI_UTILSKERNEL_LIBRARY scai_utilskernel
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_UTILSKERNEL_LIBRARY_PATH}
    ${SCAI_UTILSKERNEL_ROOT}/lib
)

set ( SCAI_UTILSKERNEL_FOUND FALSE )
if ( SCAI_UTILSKERNEL_INCLUDE_DIR )
    if ( SCAI_UTILSKERNEL_LIBRARY)
        set ( SCAI_UTILSKERNEL_FOUND TRUE )
    endif ( SCAI_UTILSKERNEL_LIBRARY )
endif ( SCAI_UTILSKERNEL_INCLUDE_DIR)

mark_as_advanced ( SCAI_UTILSKERNEL_FOUND SCAI_UTILSKERNEL_INCLUDE_DIR SCAI_UTILSKERNEL_LIBRARY )
