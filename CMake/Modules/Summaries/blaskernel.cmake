###
 # @file Summaries/blaskernel.cmake
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
 # @brief SCAI blaskernel Configuration Summary
 # @author Eric Schricker
 # @date 11.01.2016
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "=========================================" )
message ( STATUS "Summary of SCAI blaskernel Configuration:" )
message ( STATUS "=========================================" )

include ( Summaries/Modules/Compiler )

# blaskernel (core)
heading ( "Required core:" )

heading2 ( "Internal Libraries" "REQUIRED_FOUND" )
    found_message ( "SCAI common"    "SCAI_COMMON_FOUND"    "REQUIRED" "Version ${SCAI_COMMON_VERSION}"    )
    found_message ( "SCAI logging"   "SCAI_LOGGING_FOUND"   "REQUIRED" "Version ${SCAI_LOGGING_VERSION}"   )
    found_message ( "SCAI tracing"   "SCAI_TRACING_FOUND"   "REQUIRED" "Version ${SCAI_TRACING_VERSION}"   )
    found_message ( "SCAI tasking"   "SCAI_TASKING_FOUND"   "REQUIRED" "Version ${SCAI_TASKING_VERSION}"   )
    found_message ( "SCAI kregistry" "SCAI_KREGISTRY_FOUND" "REQUIRED" "Version ${SCAI_KREGISTRY_VERSION}" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_BLAS_FOUND )
    set( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_BLAS_FOUND )

heading2 ( "External Libraries" "REQUIRED_FOUND" )

    include ( Summaries/Modules/BLAS )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_KREGISTRY_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_KREGISTRY_FOUND )

heading ( "Optional External Libraries" )
include ( Summaries/Modules/Accelerator )

include ( Summaries/Modules/Build )
                       
include ( Summaries/Modules/Configuration )
