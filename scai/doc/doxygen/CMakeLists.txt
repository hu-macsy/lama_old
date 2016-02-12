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
 # @since 1.0.0
###
### Set required CMake version
cmake_minimum_required ( VERSION 2.8 )

### Set PROJECT NAME / supported LANGUAGES
project ( sphinx_dox )

### Set Module Path
set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../../../CMake/Modules" )

include ( Package/Doxygen )
include ( VersionDefinition )
include ( SetProjectVariables )

if    ( DOXYGEN_FOUND AND BUILD_DOC )
    
    configure_file ( "${CMAKE_SOURCE_DIR}/LAMA.Doxyfile.in" "${CMAKE_CURRENT_BINARY_DIR}/LAMA.Doxyfile" )

    file ( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html )

    # The initial rm command gets rid of everything previously built by this
    # custom command.

    add_custom_command (
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
        #COMMAND rm -rf ${CMAKE_CURRENT_BINARY_DIR}
        #COMMAND mkdir ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND ${DOXYGEN_EXECUTABLE} LAMA.Doxyfile
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/LAMA.Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )

    add_custom_target (
        doxygendoc
        DEPENDS
        ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
    )

    install(
        DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
        DESTINATION ${DOC_ROOT_DIR}/system
    )
endif ( DOXYGEN_FOUND AND BUILD_DOC )
