###
 # @file doc/CMakeLists.txt
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
 # @brief configuration of sphinx
 # @author Jan Ecker
 # @date 03.11.2015
###

if ( SPHINX_FOUND AND BUILD_DOC )
    set ( ${PROJECT_NAME}_DOC_DIR "${CMAKE_SOURCE_DIR}/doc" )
    set ( SPHINX_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SCAI_DOC_TYPE}" )

    # Sphinx configuration file conf.py will be generated with correct LAMA version, copyright, etc.
    # must be in same directory as Sphinx source files
    setIntersphinxInternalVariables()
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../conf.py.in" "${CMAKE_CURRENT_BINARY_DIR}/conf.py" )
    
    # Add custom command for building the userdoc
    add_custom_command (
        OUTPUT ${SPHINX_BINARY_DIR}/index.${DOC_EXTENTSION}
        COMMAND ${Sphinx-build_EXECUTABLE}
                    -b ${SCAI_DOC_TYPE}
                    -c ${CMAKE_CURRENT_BINARY_DIR}
                    -d ${SPHINX_BINARY_DIR}/_doctree
                    ${${PROJECT_NAME}_DOC_DIR}
                    ${SPHINX_BINARY_DIR}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py
        WORKING_DIRECTORY ${${PROJECT_NAME}_DOC_DIR}
    )
    
    # Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
    add_custom_target (
        doc ALL
        DEPENDS ${SPHINX_BINARY_DIR}/index.${DOC_EXTENTSION}
    )
    
    # Add documentation to the install target
    install(
        DIRECTORY ${SPHINX_BINARY_DIR}/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_DOC_DIR}
    )
endif ( SPHINX_FOUND AND BUILD_DOC )
