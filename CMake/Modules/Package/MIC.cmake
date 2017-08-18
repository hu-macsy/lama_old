###
 # @file Package/MIC.cmake
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
 # @brief Configuration for using MIC Intel Many Integrated Core Architecture
 # @author Thomas Brandes
 # @date 05.07.2013
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

### Multiple modules might call this exernal package

scai_pragma_once ()

# detect whether a MIC device is available, then USE_MIC
#
# Output variables:
#
#  USE_MIC   cached variable regarding use of MIC
#  MIC_FOUND is true if MIC device has been found
#

if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

	execute_process ( COMMAND micinfo
					  RESULT_VARIABLE MIC_RUN_RESULT_VAR
					  OUTPUT_VARIABLE MIC_RUN_OUTPUT_VAR)

    # message ( STATUS "MIC_RUN = ${MIC_RUN_RESULT_VAR} out = ${MIC_RUN_OUTPUT_VAR}" )

	# if micinfo return with SUCCESS MIC_RUN_RESULT_VAR is 0

	if ( NOT MIC_RUN_RESULT_VAR )
		set ( MIC_FOUND TRUE )
	endif ()

endif ()

scai_build_variable ( NAME      USE_MIC
                      BOOL 
                      DEFAULT   ${MIC_FOUND}
                      DOCSTRING "Enable / Disablue use of Intel Xeon PHI Coprocessor (MIC)" )

if ( USE_MIC AND NOT MIC_FOUND )
    message ( ERROR  "USE_MIC forced but not MIC Device found" )
endif ()

scai_summary_external ( NAME      MIC
                        ENABLED   ${USE_MIC} 
                        FOUND     ${MIC_FOUND}  )
