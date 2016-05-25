###
 # @file Package/GraphPartitioning.cmake
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
 # @brief Everything needed for using Graph partitioned Distribution ( only with Metis/ParMetis yet )
 # @author Lauretta Schubert
 # @date 19.08.2015
###

### USE_GRAPHPARTITIONING              - if Graph Partitioning is enabled
### SCAI_GRAPHPARTITIONING_INCLUDE_DIR - Graph Partitioning include directory
### SCAI_GRAPHPARTITIONING_LIBRARIES   - all needed Graph Partitioning libraries
### GRAPHPARTITIONING_ENABLED          - if METIS_FOUND (and PARMETIS_FOUND) and USE_GRAPHPARTITIONING

find_package ( Metis ${SCAI_FIND_PACKAGE_FLAGS} )
if    ( METIS_FOUND )
	find_package ( ParMetis ${SCAI_FIND_PACKAGE_FLAGS} )
endif ( METIS_FOUND )

## ALLOW to switch off GRAPHPARTITIONING explicitly ( doing something linke setAndCheckCache )
include ( Functions/setAndCheckCache )
setAndCheckCache ( METIS GRAPHPARTITIONING )
set ( USE_GRAPHPARTITIONING ${USE_GRAPHPARTITIONING} CACHE BOOL "Enable / Disable use of Graphpartitioning" )

set ( GRAPHPARTITIONING_ENABLED FALSE )
if    ( METIS_FOUND )
	## get Metis version
	try_run ( METIS_RUN_RESULT_VAR METIS_COMPILE_RESULT_VAR
	    ${CMAKE_BINARY_DIR}/VersionCheck
	    ${CMAKE_MODULE_PATH}/VersionCheck/metis.cpp
	    CMAKE_FLAGS 
	    -DINCLUDE_DIRECTORIES:STRING=${METIS_INCLUDE_DIR}
	    COMPILE_OUTPUT_VARIABLE METIS_COMPILE_OUTPUT_VAR
	    RUN_OUTPUT_VARIABLE METIS_RUN_OUTPUT_VAR )

	    set ( METIS_VERSION ${METIS_RUN_OUTPUT_VAR} )
	## end get Metis version

	set ( SCAI_GRAPHPARTITIONING_INCLUDE_DIR ${METIS_INCLUDE_DIR} )
	set ( SCAI_GRAPHPARTITIONING_LIBRARIES ${METIS_LIBRARY} )
	if    ( PARMETIS_FOUND AND MPI_FOUND )
		## get ParMetis version
		try_run ( PARMETIS_RUN_RESULT_VAR PARMETIS_COMPILE_RESULT_VAR
		    ${CMAKE_BINARY_DIR}/VersionCheck
		    ${CMAKE_MODULE_PATH}/VersionCheck/parmetis.cpp
		    CMAKE_FLAGS 
		    -DINCLUDE_DIRECTORIES:STRING=${PARMETIS_INCLUDE_DIR}\;${METIS_INCLUDE_DIR}\;${SCAI_MPI_INCLUDE_DIR}
		    LINK_LIBRARIES "${SCAI_MPI_LIBRARIES}"
		    COMPILE_OUTPUT_VARIABLE PARMETIS_COMPILE_OUTPUT_VAR
		    RUN_OUTPUT_VARIABLE PARMETIS_RUN_OUTPUT_VAR )

		    set ( PARMETIS_VERSION ${PARMETIS_RUN_OUTPUT_VAR} )
		## end get ParMetis version

		list ( APPEND SCAI_GRAPHPARTITIONING_INCLUDE_DIR ${PARMETIS_INCLUDE_DIR} )
		list ( APPEND SCAI_GRAPHPARTITIONING_LIBRARIES ${PARMETIS_LIBRARY} )
	endif ( PARMETIS_FOUND AND MPI_FOUND )

	if    ( USE_GRAPHPARTITIONING )
		set ( GRAPHPARTITIONING_ENABLED TRUE )
	endif ( USE_GRAPHPARTITIONING )

endif ( METIS_FOUND )