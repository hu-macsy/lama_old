###
 # @file CMake/Modules/Functions/scaiFunctions.cmake
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
 # @brief macros for adding all packages, include dirs, link libraries depending
 #        on internal and external project dependencies
 # @author Lauretta Schubert
 # @date 05.02.2016
###

macro    ( addInternalPackages )
	foreach    ( PACKAGE_TO_FIND ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
    	find_package ( ${PACKAGE_TO_FIND} ${SCAI_FIND_PACKAGE_FLAGS} REQUIRED )
	endforeach ( PACKAGE_TO_FIND ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
endmacro ( addInternalPackages )

macro    ( addExternalPackages )
	foreach    ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} Sphinx )
    	include( Package/${module} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} Sphinx )
endmacro ( addExternalPackages )

macro    ( addInternalAndExternalPackages )
	addInternalPackages()
	addExternalPackages()
endmacro ( addInternalAndExternalPackages )

## adding includes dirs from packages

macro    ( addInternalIncludes )
	foreach    ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	include_directories( ${${upper_module}_INCLUDE_DIR} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
endmacro ( addInternalIncludes )

macro    ( addExternalIncludes )
	foreach    ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	include_directories( ${SCAI_${upper_module}_INCLUDE_DIR} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
endmacro ( addExternalIncludes )

macro    ( addInternalAndExternalIncludes )
	addInternalIncludes()
	addExternalIncludes()
endmacro ( addInternalAndExternalIncludes )

## adding link libraries of packages

macro    ( addInternalLinkLibraries )
	set ( REVERT_LIST ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} ) # because list does not accept variable recursion
	list ( REVERSE REVERT_LIST )
	foreach    ( module ${REVERT_LIST} )
		string ( TOUPPER ${module} upper_module )
		set ( ${UPPER_PROJECT_NAME}_LINK_LIBRARIES ${${UPPER_PROJECT_NAME}_LINK_LIBRARIES} ${${upper_module}_LIBRARY} )
	endforeach ( module ${REVERT_LIST} )
	target_link_libraries ( ${PROJECT_NAME} ${SCAI_START_LINK_LIBRARIES} ${${UPPER_PROJECT_NAME}_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES} )
endmacro ( addInternalLinkLibraries )

macro    ( addExternalLinkLibraries )
	foreach    ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
		string ( TOUPPER ${module} upper_module )
    	target_link_libraries ( ${PROJECT_NAME} ${SCAI_${upper_module}_LIBRARIES} )
	endforeach ( module ${${UPPER_PROJECT_NAME}_EXTERNAL_DEPS} )
endmacro ( addExternalLinkLibraries )

macro    ( addInternalAndExternalLinkLibraries )
	addInternalLinkLibraries()
	addExternalLinkLibraries()
endmacro ( addInternalAndExternalLinkLibraries )

macro    ( setIntersphinxInternalVariables )
	## set PROJECT_DOC_DIR's
	foreach    ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
	    string ( TOUPPER ${module} upper_module )
	    string ( SUBSTRING ${module} 5 -1 MODULE_SURNAME )
	    set ( ${upper_module}_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}/scai-${MODULE_SURNAME}-${${upper_module}_VERSION}/ )
	endforeach ( module ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )

	## set intersphinx mapping string
	foreach    ( PACKAGE_TO_FIND ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )
		string ( TOUPPER ${PACKAGE_TO_FIND} UPPER_PACKAGE_TO_FIND )
		string ( SUBSTRING ${PACKAGE_TO_FIND} 5 -1 MODULE_SURNAME )
		if    ( MAPPING ) ##append
    		set ( MAPPING "${MAPPING},\n 'scai${MODULE_SURNAME}': ('${${UPPER_PACKAGE_TO_FIND}_DOC_DIR}', None)" )
    	else  ( MAPPING ) ##start
    		set ( MAPPING "'scai${MODULE_SURNAME}': ('${${UPPER_PACKAGE_TO_FIND}_DOC_DIR}', None)" )
    	endif ( MAPPING )
	endforeach ( PACKAGE_TO_FIND ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} )

	if    ( MAPPING )
		set ( INTERSPHINX_MAPPING "intersphinx_mapping = { ${MAPPING} }")
	endif ( MAPPING )
endmacro ( setIntersphinxInternalVariables )
