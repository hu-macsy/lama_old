###
 # @file Package/Metis.cmake
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
 # @brief Using external package METIS
 # @author Thomas Brandes
 # @date 13.07.2017
###

find_package ( Metis ${SCAI_FIND_PACKAGE_FLAGS} )

## returns METIS_FOUND

scai_build_variable ( NAME      USE_METIS
                      BOOL 
                      DEFAULT   ${METIS_FOUND}
                      DOCSTRING "use of METIS libray (graph partitioning)" )

if ( METIS_FOUND )

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

	set ( SCAI_METIS_INCLUDE_DIR ${METIS_INCLUDE_DIR} )
	set ( SCAI_METIS_LIBRARIES ${METIS_LIBRARY} )

endif ( METIS_FOUND )

scai_summary_external ( NAME      Metis
                        ENABLED   ${USE_METIS}
                        FOUND     ${METIS_FOUND} 
                        VERSION   ${METIS_VERSION} 
                        INCLUDE   ${SCAI_METIS_INCLUDE_DIR} 
                        LIBRARIES ${SCAI_METIS_LIBRARIES} )
