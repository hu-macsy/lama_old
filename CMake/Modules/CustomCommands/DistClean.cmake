###
 # @file CustomCommands/DistClean.cmake
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
 # @brief Add custom target distclean
 # @author Fraunhofer SCAI
 # @date 09.06.2015
###

if ( TARGET distclean )
    # Target already available, do no create it then anymore
else ( TARGET distclean )
    add_custom_target ( distclean )
    add_custom_command (
        TARGET distclean
        DEPENDS clean
        # make docclean (not command itself becaue it depends on clean --> doubled cmake call)
		COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/sphinx/
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/sphinx/
		# make doxygendocclean
		COMMAND ${CMAKE_COMMAND} -E remove_directory ${CMAKE_CURRENT_BINARY_DIR}/doxygen/
		COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/doxygen/
        COMMAND cd ${CMAKE_CURRENT_BINARY_DIR}
        COMMAND sh ${CMAKE_MODULE_PATH}/CustomCommands/distclean.sh
    )
endif ( TARGET distclean )
