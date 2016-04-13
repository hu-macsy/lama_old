#
# Find the kregistry (kernel registry) includes and libraries
#
# SCAI_BLASKERNEL_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_BLASKERNEL_INCLUDE_DIR - the memory include dir
# SCAI_BLASKERNEL_LIBRARY     - libraries to link against

if ( NOT SCAI_BLASKERNEL_INCLUDE_DIR )
    find_path ( SCAI_BLASKERNEL_INCLUDE_DIR blaskernel.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_BLASKERNEL_INCLUDE_PATH}/scai
        ${SCAI_BLASKERNEL_ROOT}/include/scai
    )
    get_filename_component ( SCAI_BLASKERNEL_INCLUDE_DIR ${SCAI_BLASKERNEL_INCLUDE_DIR} PATH )
endif ( NOT SCAI_BLASKERNEL_INCLUDE_DIR )

set ( SCAI_BLASKERNEL_INCLUDE_DIR ${SCAI_BLASKERNEL_INCLUDE_DIR} CACHE PATH "Path to BLASKERNEL include dir" FORCE )

find_library ( SCAI_BLASKERNEL_LIBRARY scai_blaskernel
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_BLASKERNEL_LIBRARY_PATH}
    ${SCAI_BLASKERNEL_ROOT}/lib
)

set ( SCAI_BLASKERNEL_FOUND FALSE )
if ( SCAI_BLASKERNEL_INCLUDE_DIR )
    if ( SCAI_BLASKERNEL_LIBRARY)
        set ( SCAI_BLASKERNEL_FOUND TRUE )
    endif ( SCAI_BLASKERNEL_LIBRARY )
endif ( SCAI_BLASKERNEL_INCLUDE_DIR)

mark_as_advanced ( SCAI_BLASKERNEL_FOUND SCAI_BLASKERNEL_INCLUDE_DIR SCAI_BLASKERNEL_LIBRARY )
