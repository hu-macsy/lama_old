###
 # @file Findscai_logging.cmake
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
 # @brief Find scai_logging
 # @author Lauretta Schubert
 # @date 14.08.2015
###

#
# Find the logging includes and libraries
#
# InputVariables:
# 
# CMAKE_INSTALL_PREFIX : directory is used to find the logging installation
# LOGGING_ROOT         : installation directory where the logging library is installed 
#
# SCAI_LOGGING_INCLUDE_PATH : environment variable used to find logging include directory
# SCAI_LOGGING_LIBRARY_PATH : environment variable used to find logging include directory
#
# OutputVariables:
#
# SCAI_LOGGING_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_LOGGING_INCLUDE_DIR - the logging include dir
# SCAI_LOGGING_LIBRARY     - libraries to link against
# SCAI_LOGGING_LEVEL       - level of logging, e.g.TRACE, DEBUG, INFO
# SCAI_LOGGING_FLAG        - compile flag for logging

if ( NOT SCAI_LOGGING_INCLUDE_DIR )
    find_path ( SCAI_LOGGING_INCLUDE_DIR logging.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LOGGING_INCLUDE_PATH}/scai
        ${SCAI_LOGGING_ROOT}/include/scai
    )
    get_filename_component ( SCAI_LOGGING_INCLUDE_DIR ${SCAI_LOGGING_INCLUDE_DIR} PATH )
endif ( NOT SCAI_LOGGING_INCLUDE_DIR )

set ( SCAI_LOGGING_INCLUDE_DIR ${SCAI_LOGGING_INCLUDE_DIR} CACHE PATH "Path to LOGGING include dir" FORCE )

find_library ( SCAI_LOGGING_LIBRARY scai_logging
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_LOGGING_LIBRARY_PATH}
    ${SCAI_LOGGING_ROOT}/lib
)

set ( SCAI_LOGGING_FOUND FALSE )
if ( SCAI_LOGGING_INCLUDE_DIR )
    if ( SCAI_LOGGING_LIBRARY)
        set ( SCAI_LOGGING_FOUND TRUE )
    endif ( SCAI_LOGGING_LIBRARY )
endif ( SCAI_LOGGING_INCLUDE_DIR)

mark_as_advanced ( SCAI_LOGGING_INCLUDE_DIR SCAI_LOGGING_LIBRARY )

include ( Settings/logging )
