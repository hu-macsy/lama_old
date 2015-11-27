###
 # @file CMakeLists.cmake
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief configuration of sphinx
 # @author Jan Ecker
 # @date 03.11.2015
 # @since 1.0.1
###

#find_package ( Sphinx ${SCAI_FIND_PACKAGE_FLAGS} )

### SPHINX DOCUMENTATION ###

if ( BUILD_DOC )
    ### install ###
    set ( LAMA_DOC_DIR "${CMAKE_SOURCE_DIR}/doc" )
    set ( SPHINX_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}" )
    set ( SPHINX_INSTALL_ROOT ${CMAKE_INSTALL_PREFIX}/share/doc/scai-common-${SCAI_COMMON_VERSION} )
    
    # Sphinx configuration file conf.py will be generated with correct LAMA version, copyright, etc.
    # must be in same directory as Sphinx source files
    configure_file ( "${CMAKE_SOURCE_DIR}/doc/conf.py.in" "${SPHINX_BUILD_ROOT}/conf.py" )
    
    # Add custom command for building the userdoc
    add_custom_command (
        OUTPUT ${SPHINX_BUILD_ROOT}/index.html
        COMMAND ${Sphinx-build_EXECUTABLE} -b ${SCAI_DOC_TYPE} -c ${SPHINX_BUILD_ROOT} ${LAMA_DOC_DIR} ${SPHINX_BUILD_ROOT}
        DEPENDS ${SPHINX_BUILD_ROOT}/conf.py
        WORKING_DIRECTORY ${LAMA_DOC_DIR}
    )
    
    # Add a custom target which executes the custom command (using ALL flag to add it to the default target as well!)
    add_custom_target (
        doc ALL
        DEPENDS ${SPHINX_BUILD_ROOT}/index.html
    )
    
    # Add documentation to the install target
    install(
        DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/
        DESTINATION ${SPHINX_INSTALL_ROOT}
    )
endif ( BUILD_DOC )
