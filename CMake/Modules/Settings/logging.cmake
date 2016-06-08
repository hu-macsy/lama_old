###
 # @file CMake/Modules/Settings/logging.cmake
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
 # @brief Sets default setting for logging depending on CMAKE_BUILD_TYPE.
 # @author Thomas Brandes
 # @date 06.04.2015
###

# Using the logging library requires setting of SCAI_LOG_LEVEL_xxx at compile time
# CMake variable SCAI_LOGGING_LEVEL (cache) is used for the correct setting

# LOG_CHOICES "TRACE" "DEBUG" "INFO" "WARN" "ERROR" "OFF"

## LOGGING Level
#
#  Debug   : use -DSCAI_LOG_LEVEL_DEBUG
#  Release : use -DSCAI_LOG_LEVEL_INFO
#  
#  For serious problems: -DLOG_LEVEL_TRACE
#  For benchmarks:       -DLOG_LEVEL_OFF (or -DLOG_LEVEL_FATAL, -DLOG_LEVEL_ERROR)

include ( Functions/checkValue )

if    ( NOT SCAI_LOGGING_LEVEL )
    if ( CMAKE_BUILD_TYPE STREQUAL "Release" )
        set ( DEFAULT_LOG_LEVEL "INFO" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    elseif ( CMAKE_BUILD_TYPE STREQUAL "Debug" )
        set ( DEFAULT_LOG_LEVEL "DEBUG" )
    else ( )
        set ( DEFAULT_LOG_LEVEL "TRACE" )
    endif ( )

    if    ( NOT scai_logging_FIND_QUIETLY )
        message ( "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE} implies default logging level ${DEFAULT_LOG_LEVEL}" )
    endif ( NOT scai_logging_FIND_QUIETLY )

endif ( NOT SCAI_LOGGING_LEVEL )

set ( SCAI_LOGGING_LEVEL ${DEFAULT_LOG_LEVEL} )
checkValue ( ${SCAI_LOGGING_LEVEL} "${SCAI_LOGGING_CHOICES}" )
set ( SCAI_LOGGING_LEVEL ${SCAI_LOGGING_LEVEL} CACHE STRING "Choose level of compile time logging: ${SCAI_LOGGING_CHOICES}" )

# compile flag for logging should not be put in the cache, avoids errors for wrong settings

set ( SCAI_LOGGING_FLAG "SCAI_LOG_LEVEL_${SCAI_LOGGING_LEVEL}" )
