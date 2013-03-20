#
# Find the PATOH includes and libraries
#
# PATOH_FOUND       - Do not attempt to use if "no" or undefined
# PATOH_INCLUDE_DIR - the include dir
# PATOH_LIBRARIES   - List of fully qualified libraries to link against
	
FIND_PATH(PATOH_INCLUDE_DIR patoh.h
	/usr/local/include
	/usr/include
	$ENV{PATOH_INCLUDE_PATH}
)

FIND_LIBRARY(PATOH_LIBRARIES patoh
	/usr/local/lib
	/usr/lib
	$ENV{PATOH_LIBRARY_PATH}
)
	
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( PATOH
    DEFAULT_MSG
    PATOH_INCLUDE_DIR
    PATOH_LIBRARIES
)

MARK_AS_ADVANCED( PATOH_INCLUDE_DIR PATOH_LIBRARIES )
