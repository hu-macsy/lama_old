###
 # @file prepareExamplesMakeInc.cmake
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
 # @brief Macro setting the configuration for make.inc in examples: right link libraries and -D flags
 # @author Lauretta Schubert
 # @date 02.03.2016
###

macro    ( prepareExamplesMakeInc )

	## set SCAI_EXAMPLE_LINK_LIBRARIES with own project and all dependent libraries
	set ( SCAI_EXAMPLE_LINK_LIBRARIES "-l${PROJECT_NAME}" )
	set ( REVERT_LIST ${${UPPER_PROJECT_NAME}_INTERNAL_DEPS} ) # because list does not accept variable recursion

	if ( REVERT_LIST ) # is empty for common
		list ( REVERSE REVERT_LIST )
		foreach    ( module ${REVERT_LIST} )
			set ( SCAI_EXAMPLE_LINK_LIBRARIES "${SCAI_EXAMPLE_LINK_LIBRARIES} -l${module}" )
		endforeach ( module ${REVERT_LIST} )
	endif ( REVERT_LIST )

	## set project specific SCAI_DEFINES

	if    ( SCAI_ASSERT_LEVEL )
		set ( SCAI_DEFINES "${SCAI_DEFINES} -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL}" )
	endif ( SCAI_ASSERT_LEVEL )
	
	if    ( SCAI_LOGGING_LEVEL )
		set ( SCAI_DEFINES "${SCAI_DEFINES} -DSCAI_LOG_LEVEL_${SCAI_LOGGING_LEVEL}" )
	endif ( SCAI_LOGGING_LEVEL )

	if    ( SCAI_TRACING_FLAG )
		set ( SCAI_DEFINES "${SCAI_DEFINES} -D${SCAI_TRACING_FLAG}" )
	endif ( SCAI_TRACING_FLAG )

endmacro ( prepareExamplesMakeInc )