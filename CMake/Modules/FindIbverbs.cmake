# - Find ibverbs
#
# This module looks for ibverbs support and defines the following values
#  IBVERBS_FOUND                   TRUE if IBVERBS has been found
#  IBVERBS_INCLUDE_DIR             the include path for IBVERBS
#  IBVERBS_LIBRARIES               the library to link against

find_path( IBVERBS_INCLUDE_DIR infiniband/verbs.h
    /usr/local/include
    /usr/include
    $ENV{IBVERBS_INCLUDE_PATH}
)

#message( STATUS "IBVERBS_INCLUDE_DIR: ${IBVERBS_INCLUDE_DIR}" )

FIND_LIBRARY( IBVERBS_LIBRARIES ibverbs 
    /usr/local/lib
    /usr/lib
    $ENV{IBVERBS_LIBRARY_PATH}
)

#message( STATUS "IBVERBS_LIBRARIES: ${IBVERBS_LIBRARIES}" )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( IBVERBS
    DEFAULT_MSG
    IBVERBS_INCLUDE_DIR
    IBVERBS_LIBRARIES
)
