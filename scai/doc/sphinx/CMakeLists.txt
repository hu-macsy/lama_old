###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief Main specification file for sphinx doc
 # @author Lauretta Schubert
 # @date 11.02.2016
 # @since 2.0.0
###
        
# Sphinx configuration file conf.py will be generated with correct LAMA version, copyright, etc.
# must be in same directory as Sphinx source files
setIntersphinxInternalVariables()
configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in" "${CMAKE_CURRENT_BINARY_DIR}/conf.py" )

set ( SPHINX_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/${SCAI_DOC_TYPE}" )
set ( SPHINX_DOC_THREE_DIR "${SPHINX_BINARY_DIR}/doctree" )

# Add custom command for building the userdoc
add_custom_command (
    OUTPUT ${SPHINX_BINARY_DIR}/index.${DOC_EXTENTSION}
    COMMAND ${Sphinx-build_EXECUTABLE} 
                -b ${SCAI_DOC_TYPE}
                -c ${CMAKE_CURRENT_BINARY_DIR}
                -d ${SPHINX_DOC_THREE_DIR}
                ${CMAKE_CURRENT_SOURCE_DIR}
                ${SPHINX_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/conf.py scai_solver
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
	    
# Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
add_custom_target (
    doc ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${SPHINX_BINARY_DIR} ${CMAKE_INSTALL_PREFIX}/${SPHINX_ROOT_DIR}/scai-documentation/
    DEPENDS ${SPHINX_BINARY_DIR}/index.${DOC_EXTENTSION}
    COMMENT "Creating sphinx doc."
)
