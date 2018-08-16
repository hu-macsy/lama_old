###
 # @file CustomCommands/DocClean.cmake
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
 # @brief Add custom target docclean
 # @author Fraunhofer SCAI
 # @date 09.06.2015
###

add_custom_target ( docclean )
add_custom_command (
	TARGET docclean
    DEPENDS clean
	COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/sphinx/
	COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/sphinx/
	COMMENT "Cleaning sphinx doc build dir: ${CMAKE_CURRENT_BINARY_DIR}/sphinx."
)

add_custom_target ( doxygendocclean )
add_custom_command (
	TARGET doxygendocclean
    DEPENDS clean
	COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/doxygen/
	COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/doxygen/
	COMMENT "Cleaning doxygen doc build dir: ${CMAKE_CURRENT_BINARY_DIR}/doxygen."
)