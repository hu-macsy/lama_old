###
 # @file Modules/FindIbverbs.cmake
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
 # @brief Find Ibverbs
 # @author Lauretta Schubert
 # @date 24.02.2016
###

# - Find ibverbs
#
# This module looks for ibverbs support and defines the following values
#  IBVERBS_FOUND                   TRUE if IBVERBS has been found
#  IBVERBS_INCLUDE_DIR             the include path for IBVERBS
#  IBVERBS_LIBRARIES               the library to link against

find_path( IBVERBS_INCLUDE_DIR infiniband/verbs.h
    /usr/local/include
    /usr/include
    $ENV{IBVERBS_INCLUDE_PATH}
    ${IBVERBS_ROOT}/include
)

if    ( SCAI_CMAKE_VERBOSE )
    message( STATUS "IBVERBS_INCLUDE_DIR: ${IBVERBS_INCLUDE_DIR}" )
endif ( SCAI_CMAKE_VERBOSE )

FIND_LIBRARY( IBVERBS_LIBRARIES ibverbs 
    /usr/local/lib
    /usr/lib
    $ENV{IBVERBS_LIBRARY_PATH}
    ${IBVERBS_ROOT}/lib
)

if ( SCAI_CMAKE_VERBOSE )
    message( STATUS "IBVERBS_LIBRARIES: ${IBVERBS_LIBRARIES}" )
endif ( SCAI_CMAKE_VERBOSE )

set ( IBVERBS_FOUND FALSE )

if ( IBVERBS_INCLUDE_DIR )
    if ( IBVERBS_LIBRARIES )
        set ( IBVERBS_FOUND TRUE )
    endif ( IBVERBS_LIBRARIES )
endif ( IBVERBS_INCLUDE_DIR )

mark_as_advanced( IBVERBS_INCLUDE_DIR IBVERBS_LIBRARIES )

