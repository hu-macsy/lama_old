###
 # @file Findscai_tasking.cmake
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
 # @endlicense
 #
 # @brief Find scai_tasking
 # @author Lauretta Schubert
 # @date 14.08.2015
###

#
# Find the tasking includes and libraries
#
# SCAI_TASKING_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_TASKING_INCLUDE_DIR - the tasking include dir
# SCAI_TASKING_LIBRARY     - libraries to link against

if ( NOT SCAI_TASKING_INCLUDE_DIR )
    find_path ( SCAI_TASKING_INCLUDE_DIR tasking.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_TASKING_INCLUDE_PATH}/scai
        ${SCAI_TASKING_ROOT}/include/scai
    )
    get_filename_component ( SCAI_TASKING_INCLUDE_DIR ${SCAI_TASKING_INCLUDE_DIR} PATH )
endif ( NOT SCAI_TASKING_INCLUDE_DIR )

set ( SCAI_TASKING_INCLUDE_DIR ${SCAI_TASKING_INCLUDE_DIR} CACHE PATH "Path to TASKING include dir" FORCE )

find_library ( SCAI_TASKING_LIBRARY scai_tasking
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_TASKING_LIBRARY_PATH}
    ${SCAI_TASKING_ROOT}/lib
)

set ( SCAI_TASKING_FOUND FALSE )
if ( SCAI_TASKING_INCLUDE_DIR )
    if ( SCAI_TASKING_LIBRARY)
        set ( SCAI_TASKING_FOUND TRUE )
    endif ( SCAI_TASKING_LIBRARY )
endif ( SCAI_TASKING_INCLUDE_DIR)

mark_as_advanced ( SCAI_TASKING_FOUND SCAI_TASKING_INCLUDE_DIR SCAI_TASKING_LIBRARY )
