###
 # @file FindMetis.cmake
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief Find Metis
 # @author Lauretta Schubert
 # @date 24.02.2015
###

#
# Find the METIS includes and libraries
#
# METIS is a set of serial programs for partitioning graphs, partitioning finite element meshes, 
# and producing fill reducing orderings for sparse matrices. The algorithms implemented in METIS
# are based on the multilevel recursive-bisection, multilevel k-way, and multi-constraint 
# partitioning schemes developed in the Kariypis Lab. 
#
# METIS_FOUND       - Do not attempt to use if "no" or undefined
# METIS_INCLUDE_DIR - the METIS include dir
# METIS_LIBRARIES   - List of fully qualified libraries to link against
	
if ( NOT DEFINED METIS_INCLUDE_DIR )

    find_path ( METIS_INCLUDE_DIR metis.h
    	        /usr/local/include
    	        /usr/include
    	        /usr/include/metis
    	        ${METIS_INCLUDE_PATH}
    	        $ENV{METIS_INCLUDE_PATH}
    	        ${METIS_ROOT}/include        
    	        $ENV{METIS_ROOT}/include         )
endif ()

find_library ( METIS_LIBRARY metis
	/usr/local/lib
	/usr/lib
	${METIS_LIBRARY_PATH}
	$ENV{METIS_LIBRARY_PATH}
	${METIS_ROOT}/lib
	$ENV{METIS_ROOT}/lib
)
	
if ( METIS_INCLUDE_DIR )

	add_definitions( -DHAVE_METIS_H=1 )
	
	if ( METIS_LIBRARY )
		set( METIS_LIBRARIES ${METIS_LIBRARY})
		set( METIS_FOUND TRUE )
	endif ()

endif ()
	
mark_as_advanced( METIS_INCLUDE_DIR METIS_LIBRARIES METIS_LIBRARY )
