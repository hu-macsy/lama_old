 # - Find ZOLTAN
 #
 # This module looks for ZOLTAN support and defines the following values
 #  ZOLTAN_FOUND                   TRUE if ZOLTAN has been found
 #  ZOLTAN_INCLUDE_DIR             the include path for ZOLTAN
 #  ZOLTAN_LIBRARIES                 the library to link against

FIND_PATH(ZOLTAN_INCLUDE_DIR zoltan.h
	/usr/local/include
	/usr/include
	$ENV{ZOLTAN_INCLUDE_PATH}
)
   
FIND_LIBRARY( ZOLTAN_LIBRARIES zoltan 
	/usr/local/lib
	/usr/lib
	$ENV{ZOLTAN_LIBRARY_PATH}
)
   
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( ZOLTAN
    DEFAULT_MSG
    ZOLTAN_INCLUDE_DIR
    ZOLTAN_LIBRARIES
)

MARK_AS_ADVANCED( ZOLTAN_INCLUDE_DIR ZOLTAN_LIBRARIES )
