 # - Find SCOTCH
 #
 # This module looks for SCOTCH support and defines the following values
 #  SCOTCH_FOUND                 TRUE if SCOTCH has been found
 #  SCOTCH_INCLUDE_DIR           the include path for SCOTCH
 #  SCOTCH_LIBRARIES             the library to link against
   
FIND_PATH(SCOTCH_INCLUDE_DIR scotch.h
	/usr/local/include
	/usr/include
	$ENV{SCOTCH_INCLUDE_PATH}
)
   
FIND_LIBRARY( SCOTCH_LIBRARY scotch 
	/usr/local/lib
	/usr/lib
	$ENV{SCOTCH_LIBRARY_PATH}
)
   
FIND_LIBRARY( SCOTCH_scotcherr_LIBRARY scotcherr 
	/usr/local/lib
	/usr/lib
	$ENV{SCOTCH_LIBRARY_PATH}
)
  
FIND_LIBRARY( SCOTCH_scotcherrexit_LIBRARY scotcherrexit 
	/usr/local/lib
	/usr/lib
	$ENV{SCOTCH_LIBRARY_PATH}
)

IF(SCOTCH_LIBRARY AND SCOTCH_scotcherr_LIBRARY)
	SET( SCOTCH_LIBRARIES ${SCOTCH_LIBRARY} ${SCOTCH_scotcherr_LIBRARY} )
ENDIF(SCOTCH_LIBRARY AND SCOTCH_scotcherr_LIBRARY)

IF(SCOTCH_LIBRARY AND SCOTCH_scotcherrexit_LIBRARY)
	SET( SCOTCH_LIBRARIES ${SCOTCH_LIBRARY} ${SCOTCH_scotcherrexit_LIBRARY} )
ENDIF(SCOTCH_LIBRARY AND SCOTCH_scotcherrexit_LIBRARY)

INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( SCOTCH
    DEFAULT_MSG
    SCOTCH_INCLUDE_DIR
    SCOTCH_LIBRARIES
)

MARK_AS_ADVANCED( SCOTCH_INCLUDE_DIR SCOTCH_LIBRARIES )
