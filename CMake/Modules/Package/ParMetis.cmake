###
 # @file Package/ParMetis.cmake
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Using external package PARMETIS
 # @author Thomas Brandes
 # @date 24.08.2017
###

find_package ( ParMetis ${SCAI_FIND_PACKAGE_FLAGS} )

## returns METIS_FOUND

scai_build_variable ( NAME      USE_PARMETIS
                      BOOL 
                      DEFAULT   ${PARMETIS_FOUND}
                      DOCSTRING "use of PARMETIS library (graph partitioning)" )

if ( PARMETIS_FOUND )

	## get Metis version

	try_run ( PARMETIS_RUN_RESULT_VAR PARMETIS_COMPILE_RESULT_VAR
	    ${CMAKE_BINARY_DIR}/VersionCheck
	    ${CMAKE_MODULE_PATH}/VersionCheck/metis.cpp
	    CMAKE_FLAGS 
	    -DINCLUDE_DIRECTORIES:STRING=${PARMETIS_INCLUDE_DIR}
	    COMPILE_OUTPUT_VARIABLE PARMETIS_COMPILE_OUTPUT_VAR
	    RUN_OUTPUT_VARIABLE PARMETIS_RUN_OUTPUT_VAR )

	    set ( PARMETIS_VERSION ${PARMETIS_RUN_OUTPUT_VAR} )

	## end get Metis version

	set ( SCAI_PARMETIS_INCLUDE_DIR ${PARMETIS_INCLUDE_DIR} )
	set ( SCAI_PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} )

endif ()

scai_summary_external ( NAME      ParMetis
                        ENABLED   ${USE_PARMETIS}
                        FOUND     ${PARMETIS_FOUND} 
                        VERSION   ${PARMETIS_VERSION} 
                        INCLUDE   ${SCAI_PARMETIS_INCLUDE_DIR} 
                        LIBRARIES ${SCAI_PARMETIS_LIBRARIES} )
