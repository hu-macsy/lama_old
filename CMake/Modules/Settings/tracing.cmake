###
 # @file CMake/Modules/Settings/tracing.cmake
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
 # @brief Sets default setting for tracing depending on CMAKE_BUILD_TYPE.
 # @author Thomas Brandes
 # @date 06.04.2015
###

# SCAI_TRACING
#
# If TRACING is disabled all SCAI_REGION macros in the code are
# ignored. Otherwise performance data can be collected
# where configuration is set at runtime via SCAI_TRACE.

if    ( DEFINED SCAI_TRACING )

    # translate arbitrary bool value to ON / OFF for nice summary

    if ( SCAI_TRACING )
        set ( SCAI_TRACING ON )
    else ( SCAI_TRACING )
        set ( SCAI_TRACING OFF )
    endif ( SCAI_TRACING )

else  ( DEFINED SCAI_TRACING )

    # Currently SCAI_TRACING disabled in Release mode, but has only minimal overhead
    # when SCAI_TRACE=OFF is set at runtime

    if    ( CMAKE_BUILD_TYPE STREQUAL "Release" )
	    set ( SCAI_TRACING OFF )
    else  ( CMAKE_BUILD_TYPE STREQUAL "Release" )
	    set ( SCAI_TRACING ON )
    endif ( CMAKE_BUILD_TYPE STREQUAL "Release" )

endif ( DEFINED SCAI_TRACING )

# SCAI_TRACING in CACHE, force it as value might be modified ( -DSCAI_TRACING=1 -> SCAI_TRACING=ON )

set ( SCAI_TRACING ${SCAI_TRACING} CACHE BOOL "Enable / Disable tracing of regions for performance analysis" FORCE )

if ( SCAI_TRACING )
    set ( SCAI_TRACING_FLAG "SCAI_TRACE_ON" )
else ( SCAI_TRACING )
    set ( SCAI_TRACING_FLAG "SCAI_TRACE_OFF" )
endif ( SCAI_TRACING )
