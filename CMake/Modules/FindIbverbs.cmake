###
 # @file FindGPI2.cmake
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

#message( STATUS "IBVERBS_INCLUDE_DIR: ${IBVERBS_INCLUDE_DIR}" )

FIND_LIBRARY( IBVERBS_LIBRARIES ibverbs 
    /usr/local/lib
    /usr/lib
    $ENV{IBVERBS_LIBRARY_PATH}
    ${IBVERBS_ROOT}/lib
)

#message( STATUS "IBVERBS_LIBRARIES: ${IBVERBS_LIBRARIES}" )

#include( FindPackageHandleStandardArgs )
#
#find_package_handle_standard_args( IBVERBS
#    DEFAULT_MSG
#    IBVERBS_INCLUDE_DIR
#    IBVERBS_LIBRARIES
#)

mark_as_advanced( IBVERBS_INCLUDE_DIR IBVERBS_LIBRARIES )