###
 # @file Functions/concatCacheArguments.cmake
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
 # @brief CMake function to concatenate cache argument after whitelisting for external dependencies of module to concatlist
 # @author Lauretta Schubert
 # @date 19.08.2015
###

macro    ( concatCacheArguments MODULE CONCATLIST )

	string ( TOUPPER ${MODULE} upper_module )

	set ( ${CONCATLIST} "" )
	set ( DEP_LIST ${${upper_module}_EXTERNAL_DEPS} )
	if     ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )
		set ( DEP_LIST ${DEP_LIST} Boost )
	endif  ( NOT CXX_SUPPORTS_C11 OR BUILD_TEST )

	foreach    ( package ${DEP_LIST} )
		string ( TOUPPER ${package} upper_package )
		set ( ${CONCATLIST} ${${CONCATLIST}} ${${upper_package}_ARGS} )
	endforeach ( package ${DEP_LIST} )

endmacro ( concatCacheArguments )
