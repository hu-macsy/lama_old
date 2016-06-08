###
 # @file Findscai_utilskernel.cmake
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
 # @brief Find scai_utilskernel
 # @author Eric Schricker
 # @date 18.02.2016
###

#
# Find the utilskernel includes and libraries
#
# SCAI_UTILSKERNEL_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_UTILSKERNEL_INCLUDE_DIR - the memory include dir
# SCAI_UTILSKERNEL_LIBRARY     - libraries to link against

if ( NOT SCAI_UTILSKERNEL_INCLUDE_DIR )
    find_path ( SCAI_UTILSKERNEL_INCLUDE_DIR utilskernel.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_UTILSKERNEL_INCLUDE_PATH}/scai
        ${SCAI_UTILSKERNEL_ROOT}/include/scai
    )
endif ( NOT SCAI_UTILSKERNEL_INCLUDE_DIR )

set ( SCAI_UTILSKERNEL_INCLUDE_DIR ${SCAI_UTILSKERNEL_INCLUDE_DIR} CACHE PATH "Path to UTILSKERNEL include dir" FORCE )

find_library ( SCAI_UTILSKERNEL_LIBRARY scai_utilskernel
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_UTILSKERNEL_LIBRARY_PATH}
    ${SCAI_UTILSKERNEL_ROOT}/lib
)

set ( SCAI_UTILSKERNEL_FOUND FALSE )
if ( SCAI_UTILSKERNEL_INCLUDE_DIR )
    if ( SCAI_UTILSKERNEL_LIBRARY)
        set ( SCAI_UTILSKERNEL_FOUND TRUE )
    endif ( SCAI_UTILSKERNEL_LIBRARY )
endif ( SCAI_UTILSKERNEL_INCLUDE_DIR)

mark_as_advanced ( SCAI_UTILSKERNEL_FOUND SCAI_UTILSKERNEL_INCLUDE_DIR SCAI_UTILSKERNEL_LIBRARY )
