###
 # @file listToString.cmake
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
 # @brief Function to return relative path between two pathes
 # @author Lauretta Schubert
 # @date 18.05.2016
###

## inspired by cmake bug report: https://cmake.org/Bug/view.php?id=14583

function    ( listToString separator input_list output_string_var )
	set ( _string "" )
	# Get list length
	list ( LENGTH input_list list_length )
	# If the list has 0 or 1 element, there is no need to loop over.
	if   ( list_length LESS 2 )
		set ( _string  "${input_list}" )
	else ( list_length LESS 2 )

		math ( EXPR last_element_index "${list_length} - 1" )
		foreach    ( index RANGE ${last_element_index} )

			# Get current item_value
			list ( GET input_list ${index} item_value )
			# .. and append to output string
			set ( _string  "${_string}${item_value}" )
			# Append separator if current element is NOT the last one.
			if    ( NOT index EQUAL last_element_index )
				set ( _string  "${_string}${separator}" )
			endif ( NOT index EQUAL last_element_index )

		endforeach ( index RANGE ${last_element_index} )

	endif ( list_length LESS 2 )
	set ( ${output_string_var} ${_string} PARENT_SCOPE )
endfunction ( listToString separator input_list output_string_var )