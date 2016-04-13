#
# Find the lama includes and libraries
#
# SCAI_SPARSEKERNEL_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_SPARSEKERNEL_INCLUDE_DIR - the memory include dir
# SCAI_SPARSEKERNEL_LIBRARY     - libraries to link against

if ( NOT SCAI_SPARSEKERNEL_INCLUDE_DIR )
    find_path ( SCAI_SPARSEKERNEL_INCLUDE_DIR sparsekernel.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LAMA_INCLUDE_PATH}/scai
        ${SCAI_LAMA_ROOT}/include/scai
    )
endif ( NOT SCAI_SPARSEKERNEL_INCLUDE_DIR )

set ( SCAI_SPARSEKERNEL_INCLUDE_DIR ${SCAI_SPARSEKERNEL_INCLUDE_DIR} CACHE PATH "Path to SparseKernel include dir" FORCE )

find_library ( SCAI_SPARSEKERNEL_LIBRARY scai_sparsekernel
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_SPARSEKERNEL_LIBRARY_PATH}
    ${SCAI_SPARSEKERNEL_ROOT}/lib
)

set ( SCAI_SPARSEKERNEL_FOUND FALSE )
if ( SCAI_SPARSEKERNEL_INCLUDE_DIR )
    if ( SCAI_SPARSEKERNEL_LIBRARY)
        set ( SCAI_SPARSEKERNEL_FOUND TRUE )
    endif ( SCAI_SPARSEKERNEL_LIBRARY )
endif ( SCAI_SPARSEKERNEL_INCLUDE_DIR)

mark_as_advanced ( SCAI_SPARSEKERNEL_FOUND SCAI_SPARSEKERNEL_INCLUDE_DIR SCAI_SPARSEKERNEL_LIBRARY )
