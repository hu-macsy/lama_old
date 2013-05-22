###
 # @file PackageDoxygen.cmake
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
 # @brief findPackage and configuration of doxygen
 # @author Jan Ecker
 # @date 25.04.2013
 # $Id$
###

### DOXYGEN DOCUMENTATION ###

if ( DOXYGEN_FOUND )
    ### install ###
    set ( LAMA_DOC_DIR "${LAMA_SOURCE_DIR}/doc/" )
    set ( DOXYGEN_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/doc" )
    set ( DOXYGEN_INSTALL_ROOT ${CMAKE_INSTALL_PREFIX})
    file ( MAKE_DIRECTORY ${DOXYGEN_BUILD_ROOT} )
    
    configure_file( "${CMAKE_SOURCE_DIR}/doc/LAMA.Doxyfile.in" "${CMAKE_CURRENT_BINARY_DIR}/doc/LAMA.Doxyfile" )

   # The initial rm command gets rid of everything previously built by this
   # custom command.

   add_custom_command (
      OUTPUT ${DOXYGEN_BUILD_ROOT}/html/index.html
      #COMMAND rm -rf ${DOXYGEN_BUILD_ROOT}
      #COMMAND mkdir ${DOXYGEN_BUILD_ROOT}
      COMMAND ${DOXYGEN_EXECUTABLE} LAMA.Doxyfile
      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/doc/LAMA.Doxyfile
      WORKING_DIRECTORY ${DOXYGEN_BUILD_ROOT}
   )

   add_custom_target (
      doc
      DEPENDS
      ${DOXYGEN_BUILD_ROOT}/html/index.html
   )
    add_custom_target (
        install_doc
        COMMAND mkdir -p ${DOXYGEN_INSTALL_ROOT}/share/doc/lama-${LAMA_VERSION}/system
        COMMAND cp -r ${DOXYGEN_BUILD_ROOT}/system/html/* ${DOXYGEN_INSTALL_ROOT}/share/doc/lama-${LAMA_VERSION}/system
        DEPENDS ${DOXYGEN_BUILD_ROOT}/system/html/index.html
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

else ( DOXYGEN_FOUND )
    if ( LAMA_CMAKE_VERBOSE )
        message ( STATUS "Not building system documentation because Doxygen not found." )
    endif ( LAMA_CMAKE_VERBOSE )
endif ( DOXYGEN_FOUND )