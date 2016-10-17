###
 # @file CMakeLists.txt
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
 #
 # Other Usage
 # Alternatively, this file may be used in accordance with the terms and
 # conditions contained in a signed written agreement between you and
 # Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 # @endlicense
 #
 # @brief CMake configuration file for lama - parallel I/O examples
 # @author Thomas Brandes
 # @date 03.07.2016
###

enable_language ( C )

find_package( OpenGL )
find_package( GLUT )

if    ( OPENGL_FOUND AND GLUT_FOUND )

	include_directories( ${OPENGL_INCLUDE_DIR} ${GLUT_INCLUDE_DIR} )

	set ( EXECUTABLES nbody test )

	foreach ( executable ${EXECUTABLES} )
	    add_executable( ${executable}.exe ${executable}.cpp )
	    target_link_libraries( ${executable}.exe ${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} )
	    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${executable}.cpp )
	    set ( EXAMPLE_EXECUTABLES "${EXAMPLE_EXECUTABLES} ${executable}.exe" )
	endforeach ( executable )

	# configure make.inc for examples
	include ( Functions/prepareExamplesMakeInc )
	prepareExamplesMakeInc()
	configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_make.inc.in" make.inc )
	configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/../../../examples_Makefile.in" install_Makefile COPYONLY )

	## install examples

	set ( PROJECT_EXAMPLE_SUB_DIR ${PROJECT_EXAMPLE_DIR}/mandelbrot )

	install ( FILES       ${EXAMPLE_FILES}
	          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
	install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/make.inc
	          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR} )
	install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/install_Makefile
	          DESTINATION ${PROJECT_EXAMPLE_SUB_DIR}
	          RENAME      Makefile  )

endif ( OPENGL_FOUND AND GLUT_FOUND )
