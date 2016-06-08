###
 # @file testMICfound.cmake
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
 # @brief Detect MIC availability with test program
 # @author Lauretta Schubert
 # @date 04.04.2016
###

if ( CMAKE_CXX_COMPILER_ID MATCHES Intel AND USE_OPENMP )

	execute_process ( COMMAND micinfo
					  RESULT_VARIABLE MIC_RUN_RESULT_VAR
					  OUTPUT_VARIABLE MIC_RUN_OUTPUT_VAR)

	# if micinfo return with SUCCESS MIC_RUN_RESULT_VAR is 0
	if    ( NOT MIC_RUN_RESULT_VAR )
		set ( USE_MIC TRUE )
	endif ( NOT MIC_RUN_RESULT_VAR )

else  ( CMAKE_CXX_COMPILER_ID MATCHES Intel AND USE_OPENMP )
	set ( USE_MIC FALSE )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel AND USE_OPENMP )

set ( USE_MIC ${USE_MIC} CACHE BOOL "Enable / Disable use of Intel Xeon Phi (MIC)" )
