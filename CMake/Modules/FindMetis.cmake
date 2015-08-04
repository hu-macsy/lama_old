#
# Find the METIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of
# sparse matrices. It can be found at:
#       http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# METIS_FOUND       - Do not attempt to use if "no" or undefined
# METIS_INCLUDE_DIR - the METIS include dir
# METIS_LIBRARIES   - List of fully qualified libraries to link against
	
IF( NOT DEFINED METIS_INCLUDE_DIR )
    FIND_PATH(METIS_INCLUDE_DIR metis.h
    	/usr/local/include
    	/usr/include
    	/usr/include/metis
    	$ENV{METIS_INCLUDE_PATH}
    	${METIS_ROOT}/include
    )
ENDIF( NOT DEFINED METIS_INCLUDE_DIR )

FIND_LIBRARY(METIS_LIBRARY metis
	/usr/local/lib
	/usr/lib
	$ENV{METIS_LIBRARY_PATH}
	${METIS_ROOT}/lib
)
	
IF(METIS_INCLUDE_DIR)
	ADD_DEFINITIONS( -DHAVE_METIS_H=1 )
	
	IF(METIS_LIBRARY)
		SET( METIS_LIBRARIES ${METIS_LIBRARY})
		SET( METIS_FOUND TRUE )
	ENDIF(METIS_LIBRARY)
ENDIF(METIS_INCLUDE_DIR)
	
MARK_AS_ADVANCED( METIS_FOUND METIS_INCLUDE_DIR METIS_LIBRARIES METIS_LIBRARY )