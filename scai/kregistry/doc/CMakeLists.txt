###
 # @file kregistry/doc/CMakeLists.txt
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
    # TODO: need to get rid of these extra links
    # set additional links to none dependent libs
    #set ( SCAI_TRACING_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}/scai-tracing-${SCAI_TRACING_VERSION}/ )
    #set ( SCAI_TASKING_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}/scai-tasking-${SCAI_TASKING_VERSION}/ )
    #set ( SCAI_HMEMO_DOC_DIR ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}/scai-hmemo-${SCAI_HMEMO_VERSION}/ )
    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../conf.py.in" "${CMAKE_CURRENT_BINARY_DIR}/conf.py" )
    
    set ( OUTPUT_NAME ${SPHINX_BINARY_DIR}/index.${DOC_EXTENSION} )

    # Add custom command for building the userdoc
    add_custom_command (
        OUTPUT ${OUTPUT_NAME}
        COMMAND ${Sphinx-build_EXECUTABLE}
                    -b ${SCAI_DOC_TYPE}
                    -c ${CMAKE_CURRENT_BINARY_DIR}
                    -d ${SPHINX_BINARY_DIR}/_doctree
                    ${${PROJECT_NAME}_DOC_DIR}
                    ${SPHINX_BINARY_DIR}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py
        WORKING_DIRECTORY ${${PROJECT_NAME}_DOC_DIR}
    )

    if    ( SCAI_DOC_TYPE STREQUAL latex )
        set ( PDF_OUTPUT_NAME ${SPHINX_BINARY_DIR}/SCAI${PROJECT_SURNAME}.pdf )

        message( STATUS "do ${CMAKE_MAKE_PROGRAM} in ${CMAKE_CURRENT_BINARY_DIR}/latex" )
        add_custom_command (
            OUTPUT  ${PDF_OUTPUT_NAME}
            COMMAND ${CMAKE_MAKE_PROGRAM}
            DEPENDS ${OUTPUT_NAME}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex
        )
    endif ( SCAI_DOC_TYPE STREQUAL latex )
    
    set ( OUTPUT_DEPENDENCIES ${OUTPUT_NAME} )
    if    ( SCAI_DOC_TYPE STREQUAL latex )
        set ( OUTPUT_DEPENDENCIES ${OUTPUT_NAME} ${PDF_OUTPUT_NAME} )
    endif ( SCAI_DOC_TYPE STREQUAL latex )

    # Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
    add_custom_target (
        doc ALL
        DEPENDS ${OUTPUT_DEPENDENCIES}
    )
    
    # Add documentation to the install target
    if    ( SCAI_DOC_TYPE STREQUAL latex )
        install(
            FILES ${PDF_OUTPUT_NAME}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}
        )
    else  ( SCAI_DOC_TYPE STREQUAL latex )
        install(
            DIRECTORY ${SPHINX_BINARY_DIR}/
            DESTINATION ${CMAKE_INSTALL_PREFIX}/${PROJECT_DOC_DIR}
        )
    endif ( SCAI_DOC_TYPE STREQUAL latex )

endif ( SPHINX_FOUND AND BUILD_DOC )
