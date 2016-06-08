###
 # @file Findscai_tracing.cmake
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
 # @brief Find scai_tracing
 # @author Lauretta Schubert
 # @date 14.08.2015
###

#
# Find the common includes and libraries
#
# SCAI_TRACING_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_TRACING_INCLUDE_DIR - the common include dir
# SCAI_TRACING_LIBRARY     - libraries to link against
# SCAI_TRACING             - can be used to disable tracing already at compile time
# SCAI_TRACING_FLAG        - flag to be set for compilaton of instrumented code

if ( NOT SCAI_TRACING_INCLUDE_DIR )
    find_path ( SCAI_TRACING_INCLUDE_DIR tracing.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_TRACING_INCLUDE_PATH}/scai
        ${SCAI_TRACING_ROOT}/include/scai
    )
    get_filename_component ( SCAI_TRACING_INCLUDE_DIR ${SCAI_TRACING_INCLUDE_DIR} PATH )
endif ( NOT SCAI_TRACING_INCLUDE_DIR )

set ( SCAI_TRACING_INCLUDE_DIR ${SCAI_TRACING_INCLUDE_DIR} CACHE PATH "Path to TRACING include dir" FORCE )

find_library ( SCAI_TRACING_LIBRARY scai_tracing
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_TRACING_LIBRARY_PATH}
    ${SCAI_TRACING_ROOT}/lib
)

set ( SCAI_TRACING_FOUND FALSE )
if ( SCAI_TRACING_INCLUDE_DIR )
    if ( SCAI_TRACING_LIBRARY)
        set ( SCAI_TRACING_FOUND TRUE )
    endif ( SCAI_TRACING_LIBRARY )
endif ( SCAI_TRACING_INCLUDE_DIR)

mark_as_advanced ( SCAI_TRACING_INCLUDE_DIR SCAI_TRACING_LIBRARY )

include ( Settings/tracing )
