###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief Everything needed for using Graph partitioned Distribution ( only with Metis/ParMetis yet )
 # @author Lauretta Schubert
 # @date 19.08.2015
 # @since 2.0.0
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

set ( GRAPHPARTITIONING_ENABLED FALSE )
if    ( METIS_FOUND )
	## get Metis version
	try_run ( METIS_RUN_RESULT_VAR METIS_COMPILE_RESULT_VAR
	    ${CMAKE_BINARY_DIR}
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
		    ${CMAKE_BINARY_DIR}
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