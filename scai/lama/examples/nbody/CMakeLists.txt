###
 # @file CMakeLists.txt
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
 # @brief CMake configuration file for lama - parallel I/O examples
 # @author Thomas Brandes
 # @date 03.07.2016
###

include ( scai_macro/scai_add_example )
include ( scai_macro/scai_example_makefile )
include ( scai_macro/scai_run_script )

enable_language ( C )

find_package( OpenGL )
find_package( GLUT )

if    ( OPENGL_FOUND AND GLUT_FOUND )

    include_directories( ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} )

    set ( EXECUTABLES nbody test )

    foreach ( executable ${EXECUTABLES} )

        scai_add_example ( EXECUTABLE ${executable}.exe
                           FILES      ${executable}.cpp )

        target_link_libraries( ${executable}.exe ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )

    endforeach ( executable )

    # use the recommended installation directory for this example directory
    scai_example_directory ( INSTALL_EXAMPLE_DIR )

    # generate build file in example directory
    scai_example_makefile ( DESTINATION ${INSTALL_EXAMPLE_DIR} )

    ## install examples sources
    install ( FILES       ${EXAMPLE_FILES}
              DESTINATION ${INSTALL_EXAMPLE_DIR} )


endif ( OPENGL_FOUND AND GLUT_FOUND )
