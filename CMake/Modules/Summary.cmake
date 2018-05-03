###
 # @file Summary.cmake
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
 # @brief scai Summary for build configuration
 # @author Jan Ecker
 # @date 25.04.2013
###

include ( scai_function/scaiMessages )
include ( scai_function/listToString )

emptyline()
message ( STATUS "==============================" )
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )
emptyline()

heading ( "External Software/Packages:" )

## SCAI_SUMMARY has been defined by the module projects via scai_summary

foreach ( item ${SCAI_SUMMARY} )
    message ( STATUS ${item} )
endforeach ()

emptyline()
heading ( "Configuration Details:" )
emptyline()

set ( SCAI_UNUSED_MODULES ${SCAI_ALL_MODULES} )
foreach ( module ${SCAI_DEFINED_MODULES} )
   list ( REMOVE_ITEM SCAI_UNUSED_MODULES ${module} )
endforeach()


listToString ( ", " "${SCAI_HOST_TYPES_LIST}" INST_LIST )

listToString ( ", " "${SCAI_USED_MODULES}" _List )
formatText( SET_LIST "${_List}" TextBlue )
listToString ( ", " "${SCAI_DEFINED_MODULES}" _List )
formatText( USED_LIST "${_List}" TextGreen )
listToString ( ", " "${SCAI_UNUSED_MODULES}" _List )
formatText( UNUSED_LIST "${_List}" TextRed )

indent_message ( "1" "SCAI ${MODULE_NAME} Version ${SCAI_VERSION}" )

emptyline()

indent_message ( "1" "Set SCAI modules    : ${SET_LIST}" )
indent_message ( "1" "Used SCAI modules   : ${USED_LIST}" )
indent_message ( "1" "Unused SCAI modules : ${UNUSED_LIST}" )
emptyline()
indent_message ( "1" "Build Type          : ${CMAKE_BUILD_TYPE}" )
indent_message ( "1" "Install Prefix      : ${CMAKE_INSTALL_PREFIX}" )
indent_message ( "1" "Library Type        : ${SCAI_LIBRARY_TYPE}" )
indent_message ( "1" "Numeric Types       : ${INST_LIST}" )
indent_message ( "1" "IndexType           : ${SCAI_INDEX_TYPE}" )
indent_message ( "1" "ASSERT Level        : ${SCAI_ASSERT_LEVEL} ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )" )
indent_message ( "1" "LOG Level           : ${SCAI_LOG_LEVEL} ( -DSCAI_LOG_LEVEL_${SCAI_LOG_LEVEL} )" ) 
indent_message ( "1" "TRACE               : ${SCAI_TRACE} ( -DSCAI_TRACE_${SCAI_TRACE} )" )
indent_message ( "1" "CODE COVERAGE       : ${USE_CODE_COVERAGE}" )

emptyline()

