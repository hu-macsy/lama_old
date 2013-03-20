 # - Find OSHMEM
 #
 # This module looks for openshmem support and defines the following values
 #  OSHMEM_FOUND                   TRUE if openshmem has been found
 #  OSHMEM_INCLUDE_DIR             the include path for openshmem
 #  OSHMEM_LIBRARIES               the library to link against

FIND_PATH(OSHMEM_INCLUDE_DIR mpp/shmem.h
	/usr/local/include
	/usr/include
	$ENV{OSHMEM_INCLUDE_PATH}
)
   
FIND_LIBRARY( OSHMEM_LIBRARIES openshmem 
	/usr/local/lib
	/usr/lib
	$ENV{OSHMEM_LIBRARY_PATH}
)
   
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( OSHMEM
    DEFAULT_MSG
    OSHMEM_INCLUDE_DIR
    OSHMEM_LIBRARIES
)

MARK_AS_ADVANCED( OSHMEM_INCLUDE_DIR OSHMEM_LIBRARIES )
