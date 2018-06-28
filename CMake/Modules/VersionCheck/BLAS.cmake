###
 # @file BLAS.cmake
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
 # @brief Check BLAS Version
 # @author Jan Ecker
 # @date 25.04.2013
###

### MKL

if    ( SCAI_BLAS_NAME MATCHES MKL )
	try_run ( MKL_RUN_RESULT_VAR MKL_COMPILE_RESULT_VAR
		    ${CMAKE_BINARY_DIR}/VersionCheck
	    	${CMAKE_MODULE_PATH}/VersionCheck/mkl.cpp
	    	CMAKE_FLAGS 
	    	-DINCLUDE_DIRECTORIES:STRING=${SCAI_SCAI_BLAS_INCLUDE_DIR}
	    	COMPILE_OUTPUT_VARIABLE MKL_COMPILE_OUTPUT_VAR
	    	RUN_OUTPUT_VARIABLE MKL_RUN_OUTPUT_VAR )

    set ( BLAS_VERSION ${MKL_RUN_OUTPUT_VAR} )
endif ( SCAI_BLAS_NAME MATCHES MKL )

### BLAS

if    ( SCAI_BLAS_NAME MATCHES BLAS AND BLAS_blas_LIBRARY )
	if    ( APPLE )
		execute_process ( COMMAND /usr/bin/otool -L ${BLAS_blas_LIBRARY} OUTPUT_VARIABLE _blas_output )
		string ( REGEX MATCH "current version ([0-9]+[.]*[0-9]*)" __blas_output ${_blas_output} )
	        string ( REGEX MATCH "([0-9]+[.]*[0-9]*)" BLAS_VERSION ${__blas_output} )
	elseif( UNIX )
		get_filename_component ( _LIB_PATH ${BLAS_blas_LIBRARY} PATH )
		file ( GLOB _LIBRARIES ${_LIB_PATH}/libblas.so* )
		file ( GLOB _GF_LIB ${_LIB_PATH}/libblas.so*gf )

		if    ( _GF_LIB )
			list ( REMOVE_ITEM _LIBRARIES ${_GF_LIB} )
		endif ( _GF_LIB )

		## look for lib with longest name (version number)
		set ( LEN 0 )
		foreach    ( _LIB ${_LIBRARIES} )
			string ( LENGTH ${_LIB} TMP_LEN )
			if    ( ${TMP_LEN} GREATER ${LEN} )
				set ( LEN ${TMP_LEN} )
				set ( VAL ${_LIB} )
			endif ( ${TMP_LEN} GREATER ${LEN} )
		endforeach ( _LIB ${_LIBRARIES} )
	
		if    ( DEFINED VAL )
	    	string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" BLAS_VERSION ${VAL} )
	    	if    ( "${BLAS_VERSION}" STREQUAL "" )
	    		string ( REGEX MATCH "([0-9]+\\.[0-9])" BLAS_VERSION ${VAL} )
	    		if    ( "${BLAS_VERSION}" STREQUAL "" )
	    			string ( REGEX MATCH "([0-9])" BLAS_VERSION ${VAL} )
				endif ( "${BLAS_VERSION}" STREQUAL "" ) 
			endif ( "${BLAS_VERSION}" STREQUAL "" )    	
		else  ( DEFINED VAL )
		    	set ( BLAS_VERSION "unknown BLAS version")
	        endif ( DEFINED VAL )
	endif ( )
endif ( SCAI_BLAS_NAME MATCHES BLAS AND BLAS_blas_LIBRARY )
