# - Find GPI2
#
# This module looks for GPI2 support and defines the following values
#  GPI2_FOUND                   TRUE if GPI2 has been found
#  GPI2_INCLUDE_DIR             the include path for GPI2
#  GPI2_LIBRARIES               the library to link against

find_path ( GPI2_INCLUDE_DIR GASPI.h
    /usr/local/include
    /usr/include
    $ENV{GPI2_INCLUDE_PATH}
    ${GPI2_ROOT}/include
)

# message( STATUS "GPI2_INCLUDE_DIR: ${GPI2_INCLUDE_DIR}" )

find_library ( GPI2_LIBRARIES GPI2 
    /usr/local/lib
    /usr/lib
    $ENV{GPI2_LIBRARY_PATH}
    ${GPI2_ROOT}/lib
)

# message( STATUS "GPI2_LIBRARIES: ${GPI2_LIBRARIES}" )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( GPI2
    DEFAULT_MSG
    GPI2_INCLUDE_DIR
    GPI2_LIBRARIES
)

mark_as_advanced( GPI2_INCLUDE_DIR GPI2_LIBRARIES )
