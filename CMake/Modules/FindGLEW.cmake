#
# Try to find GLEW library and include path.
# Once done this will define
#
# GLEW_FOUND
# GLEW_INCLUDE_PATH
# GLEW_LIBRARY
# 

IF (WIN32)
	FIND_PATH( GLEW_INCLUDE_PATH GL/glew.h
		$ENV{PROGRAMFILES}/GLEW/include
		${PROJECT_SOURCE_DIR}/src/nvgl/glew/include
		DOC "The directory where GL/glew.h resides")
	FIND_LIBRARY( GLEW_LIBRARY
		NAMES glew GLEW glew32 glew32s
		PATHS
		$ENV{PROGRAMFILES}/GLEW/lib
		${PROJECT_SOURCE_DIR}/src/nvgl/glew/bin
		${PROJECT_SOURCE_DIR}/src/nvgl/glew/lib
		DOC "The GLEW library")
ELSE (WIN32)
	FIND_PATH( GLEW_INCLUDE_PATH GL/glew.h
	NAMES
	glew.h
	PATHS
	/usr/include
	/usr/local/include
	/sw/include
	/opt/local/include
	~/local_root/*/include
	DOC "The directory where GL/glew.h resides")
	
	FIND_LIBRARY( GLEW_LIBRARY GLEW
		NAMES
		libGLEW.so
		PATHS
		/usr/lib64
		/usr/lib
		/usr/local/lib64
		/usr/local/lib
		/sw/lib
		/opt/local/lib
		~/local_root/*/lib
		DOC "The GLEW library")
ENDIF (WIN32)

IF (GLEW_INCLUDE_PATH)
	SET( GLEW_FOUND TRUE )
		message(STATUS "Glew-Found")
		message(STATUS "Include-Path: " ${GLEW_INCLUDE_PATH})
ELSE (GLEW_INCLUDE_PATH)
	SET( GLEW_FOUND FALSE)
		message(STATUS "Glew-NOT-Found")
ENDIF (GLEW_INCLUDE_PATH)

message(STATUS "library " ${GLEW_LIBRARY})

MARK_AS_ADVANCED( GLEW_FOUND )