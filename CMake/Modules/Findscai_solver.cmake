###
 # @file Findscai_solver.cmake
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
 # @brief Find scai_solver
 # @author Lauretta Schubert
 # @date 26.01.2016
###

#
# Find the lama includes and libraries
#
# SCAI_SOLVER_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_SOLVER_INCLUDE_DIR - the memory include dir
# SCAI_SOLVER_LIBRARY     - libraries to link against

if ( NOT SCAI_SOLVER_INCLUDE_DIR )
    find_path ( SCAI_SOLVER_INCLUDE_DIR solver.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_SOLVER_INCLUDE_PATH}/scai
        ${SCAI_SOLVER_ROOT}/include/scai
    )
endif ( NOT SCAI_SOLVER_INCLUDE_DIR )

set ( SCAI_SOLVER_INCLUDE_DIR ${SCAI_SOLVER_INCLUDE_DIR} CACHE PATH "Path to SOLVER include dir" FORCE )

find_library ( SCAI_SOLVER_LIBRARY scai_solver
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_SOLVER_LIBRARY_PATH}
    ${SCAI_SOLVER_ROOT}/lib
)

set ( SCAI_SOLVER_FOUND FALSE )
if ( SCAI_SOLVER_INCLUDE_DIR )
    if ( SCAI_SOLVER_LIBRARY)
        set ( SCAI_SOLVER_FOUND TRUE )
    endif ( SCAI_SOLVER_LIBRARY )
endif ( SCAI_SOLVER_INCLUDE_DIR)

mark_as_advanced ( SCAI_SOLVER_FOUND SCAI_SOLVER_INCLUDE_DIR SCAI_SOLVER_LIBRARY )

