 # - Find GPI
 #
 # This module looks for GPI support and defines the following values
 #  GPI_FOUND                   TRUE if GPI has been found
 #  GPI_INCLUDE_DIR             the include path for GPI
 #  GPI_LIBRARIES               the library to link against

FIND_PATH(GPI_INCLUDE_DIR GPI.h
	/usr/local/include
	/usr/include
	$ENV{GPI_INCLUDE_PATH}
)
message( STATUS "GPI_INCLUDE_DIR: ${GPI_INCLUDE_DIR}" )
   
FIND_LIBRARY( GPI_LIBRARIES GPI 
	/usr/local/lib
	/usr/lib
	$ENV{GPI_LIBRARY_PATH}
)

message( STATUS "GPI_LIBRARIES: ${GPI_LIBRARIES}" )
   
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GPI
    DEFAULT_MSG
    GPI_INCLUDE_DIR
    GPI_LIBRARIES
)

MARK_AS_ADVANCED( GPI_INCLUDE_DIR GPI_LIBRARIES )
