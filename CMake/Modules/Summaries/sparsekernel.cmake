###
 # @file Summaries/sparsekernel.cmake
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
 # @brief SCAI sparsekernel Configuration Summary
 # @author Eric Schricker
 # @date 18.02.2016
###

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "===========================================" )
message ( STATUS "Summary of SCAI sparsekernel Configuration:" )
message ( STATUS "===========================================" )

include ( Summaries/Modules/Compiler )

# SparseKernel (core)
heading ( "Required core" )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_HMEMO_FOUND 
            AND SCAI_KREGISTRY_FOUND AND SCAI_UTILSKERNEL_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_HMEMO_FOUND 
            AND SCAI_KREGISTRY_FOUND AND SCAI_UTILSKERNEL_FOUND )

heading2 ( "Internal Libraries" "REQUIRED_FOUND" )
    found_message ( "SCAI common"       "SCAI_COMMON_FOUND"       "REQUIRED" "Version ${SCAI_COMMON_VERSION}"       )
    found_message ( "SCAI logging"      "SCAI_LOGGING_FOUND"      "REQUIRED" "Version ${SCAI_LOGGING_VERSION}"      )
    found_message ( "SCAI tracing"      "SCAI_TRACING_FOUND"      "REQUIRED" "Version ${SCAI_TRACING_VERSION}"      )
    found_message ( "SCAI tasking"      "SCAI_TASKING_FOUND"      "REQUIRED" "Version ${SCAI_TASKING_VERSION}"      )
    found_message ( "SCAI hmemo"        "SCAI_HMEMO_FOUND"        "REQUIRED" "Version ${SCAI_HMEMO_VERSION}"        )
    found_message ( "SCAI kregistry"    "SCAI_KREGISTRY_FOUND"    "REQUIRED" "Version ${SCAI_KREGISTRY_VERSION}"    )
    found_message ( "SCAI utilskernel"  "SCAI_UTILSKERNEL_FOUND"  "REQUIRED" "Version ${SCAI_UTILSKERNEL_VERSION}"  )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_BLAS_FOUND )
    set( REQUIRED_FOUND TRUE )
    if ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
        set( REQUIRED_FOUND FALSE )
    endif ( SCAI_BLAS_NAME MATCHES "BLAS" AND NOT LAPACK_FOUND )
endif ( SCAI_BLAS_FOUND )

heading2 ( "External Libraries" "REQUIRED_FOUND" )

    include ( Summaries/Modules/BLAS )

heading ( "Optional External Libraries" "" )
include ( Summaries/Modules/Accelerator )

set ( REQUIRED_FOUND FALSE )
if    ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_KREGISTRY_FOUND AND SCAI_HMEMO_FOUND AND SCAI_UTILSKERNEL_FOUND )
  set ( REQUIRED_FOUND TRUE )
endif ( SCAI_COMMON_FOUND AND SCAI_LOGGING_FOUND AND SCAI_TRACING_FOUND AND SCAI_TASKING_FOUND AND SCAI_KREGISTRY_FOUND AND SCAI_HMEMO_FOUND AND SCAI_UTILSKERNEL_FOUND )

include ( Summaries/Modules/Build )  

include ( Summaries/Modules/Configuration )
