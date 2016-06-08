###
 # @file CMakeLists.txt
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
 # @brief Main specification file for sphinx doc
 # @author Lauretta Schubert
 # @date 11.02.2016
###
        
# Sphinx configuration file conf.py will be generated with correct LAMA version, copyright, etc.
# must be in same directory as Sphinx source files
setIntersphinxInternalVariables()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in" "${CMAKE_CURRENT_BINARY_DIR}/conf.py" )

set ( SPHINX_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SCAI_DOC_TYPE}" )
set ( SPHINX_DOC_THREE_DIR "${SPHINX_BINARY_DIR}/doctree" )

set ( OUTPUT_NAME ${SPHINX_BINARY_DIR}/index.${DOC_EXTENSION} )

# Add custom command for building the userdoc
add_custom_command (
    OUTPUT ${OUTPUT_NAME}
    COMMAND ${Sphinx-build_EXECUTABLE} 
                -b ${SCAI_DOC_TYPE}
                -c ${CMAKE_CURRENT_BINARY_DIR}
                -d ${SPHINX_DOC_THREE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}
                ${SPHINX_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py scai_solver
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
	
if    ( SCAI_DOC_TYPE STREQUAL latex )
    set ( PDF_OUTPUT_NAME ${SPHINX_BINARY_DIR}/LAMA.pdf )

    message( STATUS "do ${CMAKE_MAKE_PROGRAM} in ${CMAKE_CURRENT_BINARY_DIR}/latex" )
    add_custom_command (
        OUTPUT  ${PDF_OUTPUT_NAME}
        COMMAND ${CMAKE_MAKE_PROGRAM}
        DEPENDS ${OUTPUT_NAME}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex
    )
endif ( SCAI_DOC_TYPE STREQUAL latex )
    
# Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
if    ( SCAI_DOC_TYPE STREQUAL latex )
    add_custom_target (
        doc ALL
        COMMAND ${CMAKE_COMMAND} -E copy ${SPHINX_BINARY_DIR}/LAMA.pdf ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}
        DEPENDS ${PDF_OUTPUT_NAME}
        COMMENT "Creating sphinx doc."
    )
else  ( SCAI_DOC_TYPE STREQUAL latex )
    add_custom_target (
        doc ALL
        COMMAND ${CMAKE_COMMAND} -E copy_directory ${SPHINX_BINARY_DIR} ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}/scai-documentation/
        DEPENDS ${OUTPUT_NAME}
        COMMENT "Creating sphinx doc."
    )
endif ( SCAI_DOC_TYPE STREQUAL latex )
