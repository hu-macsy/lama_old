#
# Find the kregistry (kernel registry) includes and libraries
#
# SCAI_KREGISTRY_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_KREGISTRY_INCLUDE_DIR - the memory include dir
# SCAI_KREGISTRY_LIBRARY     - libraries to link against

if ( NOT SCAI_KREGISTRY_INCLUDE_DIR )
    find_path ( SCAI_KREGISTRY_INCLUDE_DIR kregistry.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_KREGISTRY_INCLUDE_PATH}/scai
        ${SCAI_KREGISTRY_ROOT}/include/scai
    )
    get_filename_component ( SCAI_KREGISTRY_INCLUDE_DIR ${SCAI_KREGISTRY_INCLUDE_DIR} PATH )
endif ( NOT SCAI_KREGISTRY_INCLUDE_DIR )

set ( SCAI_KREGISTRY_INCLUDE_DIR ${SCAI_KREGISTRY_INCLUDE_DIR} CACHE PATH "Path to KREGISTRY include dir" FORCE )

find_library ( SCAI_KREGISTRY_LIBRARY scai_kregistry
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_KREGISTRY_LIBRARY_PATH}
    ${SCAI_KREGISTRY_ROOT}/lib
)

set ( SCAI_KREGISTRY_FOUND FALSE )
if ( SCAI_KREGISTRY_INCLUDE_DIR )
    if ( SCAI_KREGISTRY_LIBRARY)
        set ( SCAI_KREGISTRY_FOUND TRUE )
    endif ( SCAI_KREGISTRY_LIBRARY )
endif ( SCAI_KREGISTRY_INCLUDE_DIR)

mark_as_advanced ( SCAI_KREGISTRY_FOUND SCAI_KREGISTRY_INCLUDE_DIR SCAI_KREGISTRY_LIBRARY )
