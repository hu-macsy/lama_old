 # - Find GASNET
 #
 # This module looks for GASNET support and defines the following values
 #  GASNET_FOUND                   TRUE if GASNET has been found
 #  GASNET_LIBRARIES               the library to link against
   
#set( GASNET_INCLUDE_DIR "/usr/local/lib" ) 
   
FIND_LIBRARY( GASNET_LIBRARIES gasnet-ibv-par
	/usr/local/lib
	/usr/lib
	$ENV{GASNET_LIBRARY_PATH}
)
   
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GASNET
    DEFAULT_MSG
    GASNET_LIBRARIES
)

MARK_AS_ADVANCED( GASNET_LIBRARIES )

