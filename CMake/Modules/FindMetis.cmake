###
 # @file FindMetis.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @brief Find Metis
 # @author Lauretta Schubert
 # @date 24.02.2015
###

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