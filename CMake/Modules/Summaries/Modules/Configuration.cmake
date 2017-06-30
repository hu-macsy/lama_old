###
 # @file CMake/Modules/Summaries/Modules/Configuration.cmake
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
 # @brief Summary concerning the configuration details (build type, library type
 #        assert level, log level, tracing, code coverage).
 # @author Lauretta Schubert
 # @date 08.04.2016
###

heading ( "Configuration Details:" )
emptyline()

set ( PROJECT_TEXT "SCAI ${PROJECT_SURNAME} Version ${SCAI_${UPPER_PROJECT_SURNAME}_VERSION}" )

if    ( ${PROJECT_NAME} MATCHES "LAMA_ALL" )
	set ( PROJECT_TEXT "SCAI ${PROJECT_SURNAME} Version ${SCAI_LAMA_ALL_VERSION} ${SCAI_VERSION_NAME}" )
endif ( ${PROJECT_NAME} MATCHES "LAMA_ALL" )

include ( Functions/listToString )
listToString ( ", " "${SCAI_HOST_TYPES_LIST}" INST_LIST )

indent_message ( "1" "${PROJECT_TEXT}" )
emptyline()
indent_message ( "1" "Build Type          : ${CMAKE_BUILD_TYPE}" )
indent_message ( "1" "Library Type        : ${SCAI_LIBRARY_TYPE}" )
indent_message ( "1" "Numeric Types       : ${INST_LIST}" )
indent_message ( "1" "IndexType           : ${SCAI_INDEX_TYPE}" )
indent_message ( "1" "ASSERT Level        : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
indent_message ( "1" "LOG Level           : ${SCAI_LOGGING_LEVEL} ( -DSCAI_LOGGING_LEVEL_${SCAI_LOGGING_LEVEL} )" ) 
indent_message ( "1" "TRACING             : ${SCAI_TRACING} ( -DSCAI_TRACING_${SCAI_TRACING} )" )

if    ( USE_CODE_COVERAGE )
    indent_message ( "1" "CODE COVERAGE       : ${USE_CODE_COVERAGE}" )
endif ( USE_CODE_COVERAGE )

emptyline()

