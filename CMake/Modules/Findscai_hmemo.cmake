###
 # @file Findscai_hmemo.cmake
 #
 # @license
 # Copyright (c) 2009-2016
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the Library of Accelerated Math Applications (LAMA).
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
 # @brief Find scai_hmemo
 # @author Thomas Brandes
 # @date 18.08.2015
###

#
# Find the hmemo (hybrid memory) includes and libraries
#
# SCAI_HMEMO_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_HMEMO_INCLUDE_DIR - the memory include dir
# SCAI_HMEMO_LIBRARY     - libraries to link against

if ( NOT SCAI_HMEMO_INCLUDE_DIR )
    find_path ( SCAI_HMEMO_INCLUDE_DIR hmemo.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_HMEMO_INCLUDE_PATH}/scai
        ${SCAI_HMEMO_ROOT}/include/scai
    )
    get_filename_component ( SCAI_HMEMO_INCLUDE_DIR ${SCAI_HMEMO_INCLUDE_DIR} PATH )
endif ( NOT SCAI_HMEMO_INCLUDE_DIR )

set ( SCAI_HMEMO_INCLUDE_DIR ${SCAI_HMEMO_INCLUDE_DIR} CACHE PATH "Path to HMEMO include dir" FORCE )

find_library ( SCAI_HMEMO_LIBRARY scai_hmemo
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_HMEMO_LIBRARY_PATH}
    ${SCAI_HMEMO_ROOT}/lib
)

set ( SCAI_HMEMO_FOUND FALSE )
if ( SCAI_HMEMO_INCLUDE_DIR )
    if ( SCAI_HMEMO_LIBRARY)
        set ( SCAI_HMEMO_FOUND TRUE )
    endif ( SCAI_HMEMO_LIBRARY )
endif ( SCAI_HMEMO_INCLUDE_DIR)

mark_as_advanced ( SCAI_HMEMO_FOUND SCAI_HMEMO_INCLUDE_DIR SCAI_HMEMO_LIBRARY )
