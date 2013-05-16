###
 # @file PackageDoxygen.cmake
 #
 # @license
 # Copyright (c) 2013
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
###

### DOXYGEN DOCUMENTATION ###

find_package( "Sphinx" )

if ( SPHINX_FOUND )
    ### install ###
    set ( LAMA_DOC_DIR "${LAMA_SOURCE_DIR}/doc" )
    set ( SPHINX_BUILD_ROOT "${CMAKE_CURRENT_BINARY_DIR}/doc/user" )
    file ( MAKE_DIRECTORY ${DOXYGEN_BUILD_ROOT} )
    
    configure_file( "${CMAKE_SOURCE_DIR}/doc/LAMA.Doxyfile.in" "${CMAKE_CURRENT_BINARY_DIR}/doc/LAMA.Doxyfile" )

   # The initial rm command gets rid of everything previously built by this
   # custom command.

   add_custom_command (
      OUTPUT ${SPHINX_BUILD_ROOT}/html/index.html
      COMMAND ${Sphinx-build_EXECUTABLE} -b html -d ${SPHINX_BUILD_ROOT}/doctrees source ${SPHINX_BUILD_ROOT}/html
      WORKING_DIRECTORY ${LAMA_DOC_DIR}/user
   )

   add_custom_target (
      userdoc
      DEPENDS
      ${SPHINX_BUILD_ROOT}/html/index.html
   )

   # Install the documentation generated at "make" time.

   # install ( DIRECTORY ${DOXYGEN_BUILD_ROOT}/ DESTINATION ${DOXYGEN_BUILD_ROOT}/html )

else ( SPHINX_FOUND )
    if ( LAMA_CMAKE_VERBOSE )
        message ( STATUS "Not building user documentation because Sphinx not found." )
    endif ( LAMA_CMAKE_VERBOSE )
endif ( SPHINX_FOUND )