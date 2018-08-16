###
 # @file DoxygenDoc.cmake
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief Commands to generate custom target doxygendoc
 # @author Jan Ecker
 # @date 03.07.2017
###

if ( DOXYGEN_FOUND )

    set ( DOXYGEN_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen" )

    configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/LAMA.Doxyfile.in" "${DOXYGEN_BINARY_DIR}/LAMA.Doxyfile" )

    file ( MAKE_DIRECTORY ${DOXYGEN_BINARY_DIR}/html )

    add_custom_command (
        OUTPUT ${DOXYGEN_BINARY_DIR}/html/index.html
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_BINARY_DIR}/LAMA.Doxyfile
        DEPENDS ${DOXYGEN_BINARY_DIR}/LAMA.Doxyfile
        WORKING_DIRECTORY ${DOXYGEN_BINARY_DIR}
    )

    add_custom_target (
        doxygendoc
        DEPENDS ${DOXYGEN_BINARY_DIR}/html/index.html
        COMMENT "Creating doxygen doc."
    )

    install ( DIRECTORY   ${DOXYGEN_BINARY_DIR}/html 
              DESTINATION ${CMAKE_INSTALL_PREFIX}/share/doc/system )

else ()

    add_custom_target (
        doxygendoc
        COMMAND echo "ATTENTION: doxygen not found, cannot build system documentation" 
    )

endif ()

