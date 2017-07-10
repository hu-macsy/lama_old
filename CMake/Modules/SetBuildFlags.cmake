###
 # @file SetBuildFlags.cmake
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
 # @brief Important CMake variable definitions
 # @author Thomas Brandes, Lauretta Schubert, Jan Ecker
 # @date 16.04.2013
###

# Check if verbose mode for CMAKE is selected
if ( DEFINED SCAI_CMAKE_VERBOSE AND SCAI_CMAKE_VERBOSE )
    set ( SCAI_FIND_PACKAGE_FLAGS )
else ()
    set ( SCAI_FIND_PACKAGE_FLAGS QUIET )
endif ()

if ( SCAI_BUILD_LIB_ONLY )
    set ( BUILD_DOC OFF )
    set ( BUILD_EXAMPLES OFF )
    set ( BUILD_TEST OFF )
endif ( SCAI_BUILD_LIB_ONLY )

# Makefile outputs
set ( CMAKE_VERBOSE_MAKEFILE OFF )

## set default switches or check user input

include ( Functions/checkValue )
include ( Functions/listToString )
include ( Functions/parseBoolean )

## DOC

include( Package/doc )

# Check if doc should be build
if    ( DEFINED BUILD_DOC )
    parseBoolean( BUILD_DOC )

    if    ( BUILD_DOC AND NOT DOC_FOUND )
        message( FATAL_ERROR "Build of documentation enabled, but configuration is incomplete!")
    endif ( BUILD_DOC AND NOT DOC_FOUND )

else  ( DEFINED BUILD_DOC )
    
    if    ( DOC_FOUND )
        set ( BUILD_DOC ON )
    else  ( DOC_FOUND )
        set ( BUILD_DOC OFF )
    endif ( DOC_FOUND )

endif ( DEFINED BUILD_DOC )
set ( BUILD_DOC ${BUILD_DOC} CACHE BOOL "Enable / Disable building of doc" )

set( DOC_ENABLED OFF )
if     ( DOC_FOUND AND BUILD_DOC )
    set( DOC_ENABLED ON )
endif ( DOC_FOUND AND BUILD_DOC )

## TEST

## CODE COVERAGE


