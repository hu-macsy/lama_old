###
 # @file Package/OpenMP.cmake
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
 # @brief findPackage and configuration of OpenMP
 # @author Jan Ecker
 # @date 25.04.2013
###

### This module returns the following variables regarding OpenMP:
### OPENMP_FOUND           - if OpenMP is found
### USE_OPENMP             - if OpenMP is enabled
### OpenMP_CXX_FLAGS       - flags to be used for compiling/linking C++ code with OpenMP pragmas
### SCAI_OMP_SCHEDULE_FLAG - needed OpenMP scheduling flag
### OPENMP_VERSION         - version id

if    ( CMAKE_VERSION VERSION_LESS 2.8.7 )
	enable_language ( C )
endif ( CMAKE_VERSION VERSION_LESS 2.8.7 ) 

find_package ( OpenMP ${SCAI_FIND_PACKAGE_FLAGS} ) # sets OPENMP_FOUND, OpenMP_CXX_FLAGS

include ( Functions/setAndCheckCache )
setAndCheckCache ( OPENMP ) # sets USE_OPENMP
set ( USE_OPENMP ${USE_OPENMP} CACHE BOOL "Enable / Disable use of OpenMP" )

# LAMA irrelevant entries will be removed from cmake GUI completely
set ( OpenMP_C_FLAGS "${OpenMP_C_FLAGS}" CACHE INTERNAL "" )

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
endif ( OPENMP_FOUND )

if    ( OPENMP_FOUND AND USE_OPENMP )

	if    ( NOT SCAI_OMP_SCHEDULE )
    	set ( SCAI_OMP_SCHEDULE "static" )
	endif ( NOT SCAI_OMP_SCHEDULE )

	#### Compile/Link flag for OpenMP will be set for all source files and all targets

	set ( SCAI_OMP_SCHEDULE_FLAG "SCAI_OMP_SCHEDULE=${SCAI_OMP_SCHEDULE}" )
	
	# Note: files using omp scheduling should be compiled with the corresponding flag
	# add_definitions ( -D${SCAI_OMP_SCHEDULE_FLAG} )

endif ( OPENMP_FOUND AND USE_OPENMP )

if    ( USE_OPENMP AND NOT OPENMP_FOUND )
    message( FATAL_ERROR "Build of LAMA with OpenMP support enabled, but configuration is incomplete!")
endif ( USE_OPENMP AND NOT OPENMP_FOUND )

set ( ADDITIONAL_CXX_FLAGS_OPENMP "${OpenMP_CXX_FLAGS}" CACHE STRING "OpenMP flag (only if enabled)" )
set ( ADDITIONAL_CXX_FLAGS_NO_OPENMP "-Wno-unknown-pragmas" ) # Supress unknown pragma warnings if OpenMP is disabled
mark_as_advanced ( ADDITIONAL_CXX_FLAGS_OPENMP ADDITIONAL_CXX_FLAGS_NO_OPENMP )
