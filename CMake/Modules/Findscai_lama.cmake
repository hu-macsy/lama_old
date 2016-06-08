###
 # @file Findscai_lama.cmake
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
 # @brief Find scai_lama
 # @author Lauretta Schubert
 # @date 26.01.2016
###

#
# Find the lama includes and libraries
#
# SCAI_LAMA_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_LAMA_INCLUDE_DIR - the memory include dir
# SCAI_LAMA_LIBRARY     - libraries to link against

if ( NOT SCAI_LAMA_INCLUDE_DIR )
    find_path ( SCAI_LAMA_INCLUDE_DIR lama.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LAMA_INCLUDE_PATH}/scai
        ${SCAI_LAMA_ROOT}/include/scai
    )
endif ( NOT SCAI_LAMA_INCLUDE_DIR )

set ( SCAI_LAMA_INCLUDE_DIR ${SCAI_LAMA_INCLUDE_DIR} CACHE PATH "Path to LAMA include dir" FORCE )

find_library ( SCAI_LAMA_LIBRARY scai_lama
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_LAMA_LIBRARY_PATH}
    ${SCAI_LAMA_ROOT}/lib
)

set ( SCAI_LAMA_FOUND FALSE )
if ( SCAI_LAMA_INCLUDE_DIR )
    if ( SCAI_LAMA_LIBRARY)
        set ( SCAI_LAMA_FOUND TRUE )
    endif ( SCAI_LAMA_LIBRARY )
endif ( SCAI_LAMA_INCLUDE_DIR)

mark_as_advanced ( SCAI_LAMA_FOUND SCAI_LAMA_INCLUDE_DIR SCAI_LAMA_LIBRARY )
