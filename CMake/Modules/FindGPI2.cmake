###
 # @file FindGPI2.cmake
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
 # @brief Find GPI2
 # @author Lauretta Schubert
 # @date 24.02.2016
###

# - Find GPI2
#
# This module looks for GPI2 support and defines the following values
#  GPI2_FOUND                   TRUE if GPI2 has been found
#  GPI2_INCLUDE_DIR             the include path for GPI2
#  GPI2_LIBRARIES               the library to link against

find_path ( GPI2_INCLUDE_DIR GASPI.h
    /usr/local/include
    /usr/include
    $ENV{GPI2_INCLUDE_PATH}
    ${GPI2_ROOT}/include
)

# message( STATUS "GPI2_INCLUDE_DIR: ${GPI2_INCLUDE_DIR}" )

find_library ( GPI2_LIBRARIES GPI2 
    /usr/local/lib
    /usr/lib
    $ENV{GPI2_LIBRARY_PATH}
    ${GPI2_ROOT}/lib
)

# message( STATUS "GPI2_LIBRARIES: ${GPI2_LIBRARIES}" )

#include( FindPackageHandleStandardArgs )
#
#find_package_handle_standard_args( GPI2
#    DEFAULT_MSG
#    GPI2_INCLUDE_DIR
#    GPI2_LIBRARIES
#)

mark_as_advanced( GPI2_INCLUDE_DIR GPI2_LIBRARIES )
