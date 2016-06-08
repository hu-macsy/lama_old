###
 # @file Functions/CheckValueAtCompleteBuild.cmake
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
 # @brief CMake function for option checking if SCAI_COMPLETE_BUILD 
 # @author Lauretta Schubert
 # @date 28.09.2015
###
include ( Functions/checkValue )
include ( Settings/switchChoices )

function ( checkValueAtCompleteBuild LIBRARY )

string ( TOUPPER ${LIBRARY} UPPER_LIBRARY )
set ( VAR_NAME "SCAI_${UPPER_LIBRARY}_EXTERNAL_DEPS" )

## from SetBuildFlags

# CMAKE_BUILD_TYPE
checkValue ( ${CMAKE_BUILD_TYPE} "${CMAKE_BUILD_TYPE_CHOICES}" )

# SCAI_LIBRARY_TYPE
checkValue ( ${SCAI_LIBRARY_TYPE} "${SCAI_LIBRARY_TYPE_CHOICES}" )

## from SCAIAssert
checkValue ( ${SCAI_ASSERT_LEVEL} "${SCAI_ASSERT_CHOICES}" )

## from SCAI_BLAS
list ( FIND ${VAR_NAME} SCAI_BLAS test_var )
#message ( STATUS "test_var SCAI_BLAS ${test_var}" )
if    ( ${test_var} GREATER -1 )
	checkValue( ${SCAI_BLAS_LIBRARY} "${SCAI_BLAS_LIBRARY_CHOICES}" )
endif ( ${test_var} GREATER -1 )

if    ( CUDA_FOUND AND USE_CUDA )
	## from SetNVCCFlags
	list ( FIND ${VAR_NAME} CUDA test_var )
	#message ( STATUS "test_var CUDA ${test_var}" )
	if    ( ${test_var} GREATER -1 )
		list ( APPEND CC_CHOICES "not-found" "13" "20" "21" "30" "32" "35" "50" )
		message ( STATUS "CUDA_COMPUTE_CAPABILITY ${CUDA_COMPUTE_CAPABILITY}" )
		message ( STATUS "CC_CHOICES ${CC_CHOICES}" )
		checkValue( ${CUDA_COMPUTE_CAPABILITY} "${CC_CHOICES}" )
	endif ( ${test_var} GREATER -1 )
endif ( CUDA_FOUND AND USE_CUDA )

endfunction ( checkValueAtCompleteBuild LIBRARY )
