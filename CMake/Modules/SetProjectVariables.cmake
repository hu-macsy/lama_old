###
 # @file Variables.cmake
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
 # @brief set project name variable for flexible use in the subprojects
 # @author Lauretta Schubert
 # @date 04.02.2016
###

## upper case project name
string ( TOUPPER ${PROJECT_NAME} UPPER_PROJECT_NAME )

## get project "sur"name (everthing after scai_)
string ( SUBSTRING ${PROJECT_NAME} 5 -1 PROJECT_SURNAME )
string ( TOUPPER ${PROJECT_SURNAME} UPPER_PROJECT_SURNAME )

## upper case project name for project dependent variables
set ( ${UPPER_PROJECT_NAME}_INCLUDE_DIR include/scai/${PROJECT_SURNAME} )

## example dir
set ( PROJECT_EXAMPLE_DIR "share/examples/scai-${PROJECT_SURNAME}-${${UPPER_PROJECT_NAME}_VERSION}" )

## doc dir
set ( DOC_ROOT_DIR share/doc )
set ( SPHINX_ROOT_DIR "${DOC_ROOT_DIR}/user/${SCAI_DOC_TYPE}/scai-${SCAI_LAMA_ALL_VERSION}")
set ( PROJECT_DOC_DIR "${SPHINX_ROOT_DIR}/scai-${PROJECT_SURNAME}-${${UPPER_PROJECT_NAME}_VERSION}" )
