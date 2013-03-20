 # - Find PTSCOTCH
 #
 # This module looks for PTSCOTCH support and defines the following values
 #  PTSCOTCH_FOUND                   TRUE if PTSCOTCH has been found
 #  PTSCOTCH_INCLUDE_DIR             the include path for PTSCOTCH
 #  PTSCOTCH_LIBRARIES               the library to link against
   
FIND_PATH(PTSCOTCH_INCLUDE_DIR ptscotch.h
	/usr/local/include
	/usr/include
	/usr/include/metis
	$ENV{PTSCOTCH_INCLUDE_PATH}
)
   
FIND_LIBRARY( PTSCOTCH_LIBRARY ptscotch 
	/usr/local/lib
	/usr/lib
	$ENV{PTSCOTCH_LIBRARY_PATH}
)
   
FIND_LIBRARY( PTSCOTCH_ptscotcherr_LIBRARY ptscotcherr 
	/usr/local/lib
	/usr/lib
	$ENV{PTSCOTCH_LIBRARY_PATH}
)
  
FIND_LIBRARY( PTSCOTCH_ptscotcherrexit_LIBRARY ptscotcherrexit 
	/usr/local/lib
	/usr/lib
	$ENV{PTSCOTCH_LIBRARY_PATH}
)

IF(PTSCOTCH_LIBRARY AND PTSCOTCH_ptscotcherr_LIBRARY)
	SET( PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${PTSCOTCH_ptscotcherr_LIBRARY} )
ELSEIF(PTSCOTCH_LIBRARY AND PTSCOTCH_ptscotcherrexit_LIBRARY)
	SET( PTSCOTCH_LIBRARIES ${PTSCOTCH_LIBRARY} ${PTSCOTCH_ptscotcherrexit_LIBRARY} )
ENDIF(PTSCOTCH_LIBRARY AND PTSCOTCH_ptscotcherr_LIBRARY)

INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( PTSCOTCH
    DEFAULT_MSG
    PTSCOTCH_INCLUDE_DIR
    PTSCOTCH_LIBRARIES
)

MESSAGE( STATUS "ptscotch lib ${PTSCOTCH_LIBRARIES} inc ${PTSCOTCH_INCLUDE_DIR}" )

MARK_AS_ADVANCED( PTSCOTCH_INCLUDE_DIR PTSCOTCH_LIBRARIES )
