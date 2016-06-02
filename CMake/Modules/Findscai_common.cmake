###
 # @file Findscai_common.cmake
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
 # @brief Find scai_common
 # @author Lauretta Schubert
 # @date 14.08.2015
###

#
# Find the common includes and libraries
#
# SCAI_COMMON_FOUND       - Do not attempt to use if "no" or undefined
# SCAI_COMMON_INCLUDE_DIR - the common include dir
# SCAI_COMMON_LIBRARY     - libraries to link against
# SCAI_COMMON_FLAGS       - Compile Flags needed to be used with libcommon

if ( NOT SCAI_COMMON_INCLUDE_DIR )
    find_path ( SCAI_COMMON_INCLUDE_DIR common.hpp
        /usr/local/include/scai
        /usr/include/scai
        ${CMAKE_INSTALL_PREFIX}/include/scai
        $ENV{SCAI_COMMON_INCLUDE_PATH}/scai
        ${SCAI_COMMON_ROOT}/include/scai
    )
    get_filename_component ( SCAI_COMMON_INCLUDE_DIR ${SCAI_COMMON_INCLUDE_DIR} PATH )
endif ( NOT SCAI_COMMON_INCLUDE_DIR )

set ( SCAI_COMMON_INCLUDE_DIR ${SCAI_COMMON_INCLUDE_DIR} CACHE PATH "Path to COMMON include dir" FORCE )

find_library ( SCAI_COMMON_LIBRARY scai_common
    /usr/local/lib
    /usr/lib
    $ENV{SCAI_COMMON_LIBRARY_PATH}
    ${SCAI_COMMON_ROOT}/lib
)

set ( SCAI_COMMON_FOUND FALSE )
if ( SCAI_COMMON_INCLUDE_DIR )
    if ( SCAI_COMMON_LIBRARY)
        set ( SCAI_COMMON_FOUND TRUE )
    endif ( SCAI_COMMON_LIBRARY )
endif ( SCAI_COMMON_INCLUDE_DIR)

# set SCAI_COMMON_FLAGS for required dependencies
set ( SCAI_COMMON_FLAGS "" )
if    ( SCAI_COMMON_FOUND )

    include ( Compiler/CheckC++11 )
    set ( SCAI_COMMON_FLAGS "${SCAI_COMMON_FLAGS} ${SCAI_LANG_FLAGS}" )

    if ( NOT CXX_SUPPORTS_C11 )

        # add Boost to SCAI_COMMON_INCLUDE_DIR

        include ( Package/Boost )
        list ( APPEND SCAI_COMMON_INCLUDE_DIR ${SCAI_BOOST_INCLUDE_DIR} )

    endif ( NOT CXX_SUPPORTS_C11 )
    
    include ( Package/OpenMP )
    if    ( OPENMP_FOUND AND USE_OPENMP )
        set ( SCAI_COMMON_FLAGS "${SCAI_COMMON_FLAGS} ${OpenMP_CXX_FLAGS}" )
    endif ( OPENMP_FOUND AND USE_OPENMP )

    # remove leading and trailing whitespaces
    string ( STRIP "${SCAI_COMMON_FLAGS}" SCAI_COMMON_FLAGS )
endif ( SCAI_COMMON_FOUND)

mark_as_advanced ( SCAI_COMMON_FOUND SCAI_COMMON_INCLUDE_DIR SCAI_COMMON_LIBRARY SCAI_COMMON_FLAGS )
