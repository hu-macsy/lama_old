###
 # @file Package/OpenMP.cmake
 #
 # @license
 # Copyright (c) 2009-2018
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # This file is part of the SCAI framework LAMA.
 #
 # LAMA is free software: you can redistribute it and/or modify it under the
 # terms of the GNU Lesser General Public License as published by the Free
 # Software Foundation, either version 3 of the License, or (at your option)
 # any later version.
 #
 # LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 # WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 # FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 # more details.
 #
 # You should have received a copy of the GNU Lesser General Public License
 # along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 # @endlicense
 #
 # @brief findPackage and configuration of OpenMP
 # @author Jan Ecker
 # @date 25.04.2013
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

### This module returns the following variables regarding OpenMP:
### OPENMP_FOUND           - if OpenMP is found
### USE_OPENMP             - if OpenMP is enabled
### OpenMP_CXX_FLAGS       - flags to be used for compiling/linking C++ code with OpenMP pragmas
### OPENMP_VERSION         - version id

scai_pragma_once( OpenMP )

find_package ( OpenMP ${SCAI_FIND_PACKAGE_FLAGS} ) # sets OPENMP_FOUND, OpenMP_CXX_FLAGS

set ( OPENMP_MINIMUM_VERSION 3.0 ) # because of use of collapse

# LAMA irrelevant entries will be removed from cmake GUI completely
if ( OpenMP_C_FLAGS )
    set ( OpenMP_C_FLAGS "${OpenMP_C_FLAGS}" CACHE INTERNAL "" )
endif ( OpenMP_C_FLAGS )

## get OpenMP version
if    ( OPENMP_FOUND )
    try_run ( OPENMP_RUN_RESULT_VAR OPENMP_COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}/VersionCheck
        ${CMAKE_MODULE_PATH}/VersionCheck/openmp.cpp
        CMAKE_FLAGS 
        -DCOMPILE_DEFINITIONS:STRING=${OpenMP_CXX_FLAGS}
        COMPILE_OUTPUT_VARIABLE OPENMP_COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE OPENMP_RUN_OUTPUT_VAR )

    set ( OPENMP_VERSION ${OPENMP_RUN_OUTPUT_VAR} )

    if    ( ${OPENMP_VERSION} VERSION_LESS ${OPENMP_MINIMUM_VERSION} )
        message ( WARNING "Found OpenMP version (${OPENMP_VERSION}) of your compiler (${CMAKE_CXX_COMPILER_ID} v ${CXX_COMPILER_VERSION}) is to old - must be at least ${OPENMP_MINIMUM_VERSION}, disable OpenMP support!!!" )
        set ( OPENMP_FOUND FALSE )
    endif ()
endif ( OPENMP_FOUND )

scai_build_variable ( NAME      USE_OPENMP   
                      BOOL 
                      DEFAULT   ${OPENMP_FOUND}
                      DOCSTRING "use of OpenMP (shared memory parallelization)" )

if    ( USE_OPENMP AND NOT OPENMP_FOUND )
    message( FATAL_ERROR "Build of LAMA with OpenMP support enabled, but configuration is incomplete!")
endif ( USE_OPENMP AND NOT OPENMP_FOUND )

set ( ADDITIONAL_CXX_FLAGS_NO_OPENMP "-Wno-unknown-pragmas" CACHE STRING "ignore OpenMP pragmas (only if disabled)" )

mark_as_advanced ( ADDITIONAL_CXX_FLAGS_NO_OPENMP )

scai_summary_enabled  ( NAME      OpenMP 
                        ENABLED   ${USE_OPENMP} )

scai_summary_external ( NAME      OpenMP
                        FOUND     ${OPENMP_FOUND} 
                        VERSION   ${OPENMP_VERSION} 
                        CXX_FLAGS ${OpenMP_CXX_FLAGS} )
