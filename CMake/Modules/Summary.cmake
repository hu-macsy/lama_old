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

include ( Functions/scaiMessages )

emptyline()
message ( STATUS "==============================" )
message ( STATUS "Summary of SCAI Configuration:" )
message ( STATUS "==============================" )

heading ( "Compiler:" )

if    ( CMAKE_CXX_COMPILER AND ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR ) )
    set( REQUIRED_FOUND TRUE )
else  ( CMAKE_CXX_COMPILER AND ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR ) )
    set( REQUIRED_FOUND FALSE )
    message ( FATAL_ERROR "Compiler Configuration incomplete" )
endif ( CMAKE_CXX_COMPILER AND ( CXX_SUPPORTS_C11 OR BOOST_INCLUDE_DIR ) )

heading2 ( "Configuration" "REQUIRED_FOUND" )
    found_message ( "C++ Compiler" "CMAKE_CXX_COMPILER" "REQUIRED" "${CMAKE_CXX_COMPILER_ID} Version ${CXX_COMPILER_VERSION}" )
    found_message ( "with C++11 support" "CXX_SUPPORTS_C11" "REQUIRED" "" )

if    ( NOT CXX_SUPPORTS_C11 )
    emptyline()
    message ( STATUS "Either compiler supporting C++11 or Boost needed." )
    found_message ( "Boost" "BOOST_INCLUDE_DIR" "REQUIRED" "Version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION} at ${BOOST_INCLUDE_DIR}" )
endif ( NOT CXX_SUPPORTS_C11 )

heading ( "Supported Accelertors" )

# OpenMP usage
heading3 ( "OpenMP" "USE_OPENMP" )
    found_message ( "OpenMP" "OPENMP_VERSION" "OPTIONAL" "Version ${OPENMP_VERSION}" )
    found_message ( "compile flag" "OpenMP_CXX_FLAGS" "OPTIONAL" "${OpenMP_CXX_FLAGS}" )

# LAMA CUDA
heading3 ( "CUDA" "CUDA_ENABLED" )
    found_message ( "CUDA" "CUDA_FOUND" "OPTIONAL" "Version ${CUDA_VERSION} at ${SCAI_CUDA_INCLUDE_DIR}" )
    found_message ( "Compute Capability" "CUDA_COMPUTE_CAPABILITY" "OPTIONAL" "${CUDA_COMPUTE_CAPABILITY}" )
                           
# LAMA MIC
heading3 ( "MIC" "USE_MIC" )

heading ( "External Dependencies:" )

message ( STATUS " " )

## SCAI_SUMMARY has been defined by the module projects via scai_summary

foreach ( item ${SCAI_SUMMARY} )
    message ( STATUS ${item} )
endforeach ()

heading ( "Build options:" "" )

# EXAMPLES
heading3 ( "Examples" "BUILD_EXAMPLES" )

# TEST
heading3 ( "Test" "BOOST_TEST_ENABLED" )
found_message ( "Boost Unit Test" "Boost_UNIT_TEST_FRAMEWORK_FOUND" "OPTIONAL" "Version ${BOOST_VERSION} at ${BOOST_INCLUDE_DIR}" )

# DOC
heading3 ( "Documentation" "DOC_ENABLED" )
    found_message ( "Sphinx" "SPHINX_FOUND" "OPTIONAL" "Version ${Sphinx_VERSION_STRING} with ${Sphinx-build_EXECUTABLE}" )
    found_message ( "Doxygen" "DOXYGEN_FOUND" "OPTIONAL" "Version ${DOXYGEN_VERSION} with ${DOXYGEN_EXECUTABLE}" )


heading ( "Configuration Details:" )
emptyline()

set ( PROJECT_TEXT "SCAI ${MODULE_NAME} Version ${SCAI_VERSION}" )

set ( SCAI_UNUSED_MODULES ${SCAI_ALL_MODULES} )
foreach ( module ${SCAI_USED_MODULES} )
   list ( REMOVE_ITEM SCAI_UNUSED_MODULES ${module} )
endforeach()

include ( Functions/listToString )
listToString ( ", " "${SCAI_HOST_TYPES_LIST}" INST_LIST )
listToString ( ", " "${SCAI_USED_MODULES}" AVAIL_LIST )
listToString ( ", " "${SCAI_UNUSED_MODULES}" UNAVAIL_LIST )

indent_message ( "1" "${PROJECT_TEXT}" )
emptyline()
indent_message ( "1" "Used SCAI modules   : ${TextGreen}${AVAIL_LIST}${TextColorReset}" )
indent_message ( "1" "Unused SCAI modules : ${TextRed}${UNAVAIL_LIST}${TextColorReset}" )
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

