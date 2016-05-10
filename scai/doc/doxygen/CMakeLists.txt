###
 # @file doc/doxygen/CMakeLists.txt
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
 # @brief configuration of doxygen
 # @author Jan Ecker
 # @date 25.04.2013
###

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/LAMA.Doxyfile.in" "${CMAKE_CURRENT_BINARY_DIR}/LAMA.Doxyfile" )

file ( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html )

# The initial rm command gets rid of everything previously built by this
# custom command.

add_custom_command (
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    #COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}
    #COMMAND mkdir ${CMAKE_CURRENT_BINARY_DIR}
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/LAMA.Doxyfile
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/LAMA.Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
)

add_custom_target (
    doxygendoc
    COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/system
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/html/ ${CMAKE_INSTALL_PREFIX}/${DOC_ROOT_DIR}/system
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    COMMENT "Creating doxygen doc."
)
