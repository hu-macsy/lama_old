###
 # @file SphinxDoc.cmake
 #
 # @license
 # Copyright (c) 2009-2017
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief Commands to generate custom target doc_${MODULE_NAME}
 # @author Jan Ecker
 # @date 03.11.2015
###

## The following variables are used and must be defined:
## 
##   MODULE_NAME    ( current module name, e.g. lama, .... )
##   INTERNAL_DEPS  ( list of required modules )
## 
##   SCAI_DOC_TYPE  ( e.g. html or latex )

macro ( setIntersphinxMapping MODULES )

    ## set intersphinx mapping string

    foreach ( module ${MODULES} )
        if    ( MAPPING ) ##append
            set ( MAPPING "${MAPPING},\n 'scai${module}': ('${CMAKE_BINARY_DIR}/${module}/doc/html', None)" )
        else  ( MAPPING ) ##start
            set ( MAPPING "'scai${module}': ('${CMAKE_BINARY_DIR}/${module}/doc/html', None)" )
        endif ( MAPPING )
    endforeach ()

    if ( MAPPING )
        set ( INTERSPHINX_MAPPING "intersphinx_mapping = { ${MAPPING} }")
    endif ()

endmacro ()

if ( SPHINX_FOUND AND BUILD_DOC )

    # message ( STATUS "make target for sphinx doc of ${MODULE_NAME}: depends on ${INTERNAL_DEPS}" )

    set ( SPHINX_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/doc" )
    set ( SPHINX_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc/${SCAI_DOC_TYPE}" )

    # Sphinx configuration file conf.py will be generated with correct LAMA version, copyright, etc.
    # must be in same directory as Sphinx source files

    setIntersphinxMapping ( "${INTERNAL_DEPS}" )

    configure_file ( "${CMAKE_SOURCE_DIR}/conf.py.in" "${CMAKE_CURRENT_BINARY_DIR}/doc/conf.py" )
    
    if ( SCAI_DOC_TYPE STREQUAL json )
        set ( DOC_EXTENSION "fjson" )
    elseif ( SCAI_DOC_TYPE STREQUAL html )
        set ( DOC_EXTENSION "html" )
    elseif ( SCAI_DOC_TYPE STREQUAL xml )
        set ( DOC_EXTENSION "xml" )
    elseif ( SCAI_DOC_TYPE STREQUAL latex )
        set ( DOC_EXTENSION "tex" )
    else ()
        message( ERROR "illegal SCAI_DOC_TYPE=${SCAI_DOC_TYPE}" )
    endif ( )

    set ( OUTPUT_NAME ${SPHINX_BINARY_DIR}/index.${DOC_EXTENSION} )

    # Add custom command for building the userdoc
    add_custom_command (
        OUTPUT ${OUTPUT_NAME}
        COMMAND ${Sphinx-build_EXECUTABLE}
                    -b ${SCAI_DOC_TYPE}
                    -c ${CMAKE_CURRENT_BINARY_DIR}/doc
                    -d ${SPHINX_BINARY_DIR}/_doctree
                    ${SPHINX_SOURCE_DIR}
                    ${SPHINX_BINARY_DIR}
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doc/conf.py
        WORKING_DIRECTORY ${SPHINX_SOURCE_DIR}
    )

    if ( SCAI_DOC_TYPE STREQUAL latex )
        set ( PDF_OUTPUT_NAME ${SPHINX_BINARY_DIR}/SCAI${MODULE_NAME}.pdf )

        message( STATUS "do ${CMAKE_MAKE_PROGRAM} in ${CMAKE_CURRENT_BINARY_DIR}/doc/latex" )
        add_custom_command (
            OUTPUT  ${PDF_OUTPUT_NAME}
            COMMAND ${CMAKE_MAKE_PROGRAM}
            DEPENDS ${OUTPUT_NAME}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex
        )
    endif ()
    
    set ( OUTPUT_DEPENDENCIES ${OUTPUT_NAME} )
    if    ( SCAI_DOC_TYPE STREQUAL latex )
        set ( OUTPUT_DEPENDENCIES ${OUTPUT_NAME} ${PDF_OUTPUT_NAME} )
    endif ( SCAI_DOC_TYPE STREQUAL latex )

    # Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
    add_custom_target (
        doc_${MODULE_NAME} ALL
        DEPENDS ${OUTPUT_DEPENDENCIES}
    )

    # add dependencies for all other docs to guarantee intersphinx consistency

    foreach ( dep_module ${INTERNAL_DEPS} )
        # message( STATUS "Add dependency doc_${MODULE_NAME} doc_${dep_module}" )
        add_dependencies ( doc_${MODULE_NAME} doc_${dep_module} )
    endforeach ()
    
    # Add documentation to the install target
    if    ( SCAI_DOC_TYPE STREQUAL latex )
        install(
            FILES ${PDF_OUTPUT_NAME}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/user//${SCAI_DOC_TYPE}/scai-${SCAI_LAMA_ALL_VERSION}
        )
    else  ( SCAI_DOC_TYPE STREQUAL latex )
        install(
            DIRECTORY ${SPHINX_BINARY_DIR}
            DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/user/${MODULE_NAME}
        )
    endif ( SCAI_DOC_TYPE STREQUAL latex )

endif ( SPHINX_FOUND AND BUILD_DOC )
