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

### GRAPH_PART_FOUND                 - if Graph Partitioning library is found
### USE_GRAPH_PART                   - if Graph Partitioning is enabled
### GRAPHPARTITIONING_INCLUDE_DIR    - Graph Partitioning include directory
### SCAI_GRAPHPARTITIONING_LIBRARIES - all needed Graph Partitioning libraries

find_package ( Metis ${SCAI_FIND_PACKAGE_FLAGS} )
if    ( METIS_FOUND )
	find_package ( ParMetis ${SCAI_FIND_PACKAGE_FLAGS} )
endif ( METIS_FOUND )

## ALLOW to switch off GRAPH_PART explicitly ( doing something linke setAndCheckCache )
include ( Functions/setAndCheckCache )
setAndCheckCache ( METIS GRAPH_PART )

if    ( METIS_FOUND )
	set ( GRAPHPARTITIONING_INCLUDE_DIR ${METIS_INCLUDE_DIR} )
	set ( SCAI_GRAPHPARTITIONING_LIBRARIES ${METIS_LIBRARY} )
	if    ( PARMETIS_FOUND )
		list ( APPEND GRAPHPARTITIONING_INCLUDE_DIR ${PARMETIS_INCLUDE_DIR} )
		list ( APPEND SCAI_GRAPHPARTITIONING_LIBRARIES ${PARMETIS_LIBRARY} )
	endif ( PARMETIS_FOUND )
endif ( METIS_FOUND )