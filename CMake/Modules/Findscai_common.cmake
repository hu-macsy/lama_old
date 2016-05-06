#
# Find the common includes and libraries
#
# SCAI_COMMON_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_COMMON_INCLUDE_DIR - the common include dir
# SCAI_COMMON_LIBRARY     - libraries to link against
# SCAI_COMMON_FLAGS       - Compile Flags needed to be used with libcommon

if ( NOT SCAI_COMMON_INCLUDE_DIR )
    find_path ( SCAI_COMMON_INCLUDE_DIR common.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_COMMON_INCLUDE_PATH}/scai
        ${SCAI_COMMON_ROOT}/include/scai
    )
    get_filename_component ( SCAI_COMMON_INCLUDE_DIR ${SCAI_COMMON_INCLUDE_DIR} PATH )
endif ( NOT SCAI_COMMON_INCLUDE_DIR )

set ( SCAI_COMMON_INCLUDE_DIR ${SCAI_COMMON_INCLUDE_DIR} CACHE PATH "Path to COMMON include dir" FORCE )

find_library ( SCAI_COMMON_LIBRARY scai_common
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_COMMON_LIBRARY_PATH}
    ${SCAI_COMMON_ROOT}/lib
)

set ( SCAI_COMMON_FOUND FALSE )
if ( SCAI_COMMON_INCLUDE_DIR )
    if ( SCAI_COMMON_LIBRARY)
        set ( SCAI_COMMON_FOUND TRUE )
    endif ( SCAI_COMMON_LIBRARY )
endif ( SCAI_COMMON_INCLUDE_DIR)

# set SCAI_COMMON_FLAGS for required dependencies
set ( SCAI_COMMON_FLAGS "" )
if    ( SCAI_COMMON_FOUND )

    include ( Compiler/CheckC++11 )
    set ( SCAI_COMMON_FLAGS "${SCAI_COMMON_FLAGS} ${SCAI_LANG_FLAGS}" )

    include ( Package/OpenMP )
    if    ( OPENMP_FOUND AND USE_OPENMP )
        set ( SCAI_COMMON_FLAGS "${SCAI_COMMON_FLAGS} ${OpenMP_CXX_FLAGS}" )
    endif ( OPENMP_FOUND AND USE_OPENMP )

    # remove leading and trailing whitespaces
    string ( STRIP "${SCAI_COMMON_FLAGS}" SCAI_COMMON_FLAGS )
endif ( SCAI_COMMON_FOUND)

mark_as_advanced ( SCAI_COMMON_FOUND SCAI_COMMON_INCLUDE_DIR SCAI_COMMON_LIBRARY SCAI_COMMON_FLAGS )
