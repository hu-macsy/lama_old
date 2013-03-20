#
# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of
# sparse matrices. It can be found at:
#       http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_FOUND       - Do not attempt to use if "no" or undefined
# PARMETIS_INCLUDE_DIR - the METIS include dir
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against
	
IF( NOT DEFINED PARMETIS_INCLUDE_DIR )
    FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h
    	/usr/local/include
    	/usr/include
    	/usr/include/metis
    	$ENV{PARMETIS_INCLUDE_PATH}
    	${PARMETIS_ROOT}/include
    )
ENDIF( NOT DEFINED PARMETIS_INCLUDE_DIR )

FIND_LIBRARY(PARMETIS_LIBRARY parmetis
	/usr/local/lib
	/usr/lib
	$ENV{PARMETIS_LIBRARY_PATH}
	${PARMETIS_ROOT}/lib
)
	
IF(PARMETIS_INCLUDE_DIR)
	ADD_DEFINITIONS( -DHAVE_PARMETIS_H=1 )
	
	IF(PARMETIS_LIBRARY)
		SET( PARMETIS_LIBRARIES ${PARMETIS_LIBRARY})
		SET( PARMETIS_FOUND TRUE )
	ENDIF(PARMETIS_LIBRARY)
ENDIF(PARMETIS_INCLUDE_DIR)

MARK_AS_ADVANCED( PARMETIS_FOUND PARMETIS_INCLUDE_DIR PARMETIS_LIBRARY )