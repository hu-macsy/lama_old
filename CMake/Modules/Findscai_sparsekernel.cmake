###
 # @file Findscai_sparsekernel.cmake
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
 # @brief Find scai_sparsekernel
 # @author Eric Schricker
 # @date 18.02.2016
###

#
# Find the lama includes and libraries
#
# SCAI_SPARSEKERNEL_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_SPARSEKERNEL_INCLUDE_DIR - the memory include dir
# SCAI_SPARSEKERNEL_LIBRARY     - libraries to link against

if ( NOT SCAI_SPARSEKERNEL_INCLUDE_DIR )
    find_path ( SCAI_SPARSEKERNEL_INCLUDE_DIR sparsekernel.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_LAMA_INCLUDE_PATH}/scai
        ${SCAI_LAMA_ROOT}/include/scai
    )
endif ( NOT SCAI_SPARSEKERNEL_INCLUDE_DIR )

set ( SCAI_SPARSEKERNEL_INCLUDE_DIR ${SCAI_SPARSEKERNEL_INCLUDE_DIR} CACHE PATH "Path to SparseKernel include dir" FORCE )

find_library ( SCAI_SPARSEKERNEL_LIBRARY scai_sparsekernel
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_SPARSEKERNEL_LIBRARY_PATH}
    ${SCAI_SPARSEKERNEL_ROOT}/lib
)

set ( SCAI_SPARSEKERNEL_FOUND FALSE )
if ( SCAI_SPARSEKERNEL_INCLUDE_DIR )
    if ( SCAI_SPARSEKERNEL_LIBRARY)
        set ( SCAI_SPARSEKERNEL_FOUND TRUE )
    endif ( SCAI_SPARSEKERNEL_LIBRARY )
endif ( SCAI_SPARSEKERNEL_INCLUDE_DIR)

mark_as_advanced ( SCAI_SPARSEKERNEL_FOUND SCAI_SPARSEKERNEL_INCLUDE_DIR SCAI_SPARSEKERNEL_LIBRARY )
