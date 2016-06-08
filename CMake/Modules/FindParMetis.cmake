###
 # @file FindParMetis.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Affero General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Affero General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Find ParMetis
 # @author Lauretta Schubert
 # @date 24.02.2015
###

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
