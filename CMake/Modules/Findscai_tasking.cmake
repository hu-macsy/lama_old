#
# Find the tasking includes and libraries
#
# SCAI_TASKING_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_TASKING_INCLUDE_DIR - the tasking include dir
# SCAI_TASKING_LIBRARY     - libraries to link against

if ( NOT SCAI_TASKING_INCLUDE_DIR )
    find_path ( SCAI_TASKING_INCLUDE_DIR tasking.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_TASKING_INCLUDE_PATH}/scai
        ${SCAI_TASKING_ROOT}/include/scai
    )
    get_filename_component ( SCAI_TASKING_INCLUDE_DIR ${SCAI_TASKING_INCLUDE_DIR} PATH )
endif ( NOT SCAI_TASKING_INCLUDE_DIR )

set ( SCAI_TASKING_INCLUDE_DIR ${SCAI_TASKING_INCLUDE_DIR} CACHE PATH "Path to TASKING include dir" FORCE )

find_library ( SCAI_TASKING_LIBRARY scai_tasking
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_TASKING_LIBRARY_PATH}
    ${SCAI_TASKING_ROOT}/lib
)

set ( SCAI_TASKING_FOUND FALSE )
if ( SCAI_TASKING_INCLUDE_DIR )
    if ( SCAI_TASKING_LIBRARY)
        set ( SCAI_TASKING_FOUND TRUE )
    endif ( SCAI_TASKING_LIBRARY )
endif ( SCAI_TASKING_INCLUDE_DIR)

mark_as_advanced ( SCAI_TASKING_FOUND SCAI_TASKING_INCLUDE_DIR SCAI_TASKING_LIBRARY )
