###
 # @file PackageMPI.cmake
 #
 # @license
 # Copyright (c) 2009-2013
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief findPackage and configuration of OpenMP
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

### This module returns the following variables regarding OpenMP:
###
### OPENMP_FOUND           - if OpenMP is found
### USE_OPENMP             - if OpenMP is enabled
### OpenMP_CXX_FLAGS       - flags to be used for compiling/linking C++ code with OpenMP pragmas
### SCAI_OMP_SCHEDULE_FLAG - needed OpenMP scheduling flag
### OPENMP_VERSION         - version id

if    ( CMAKE_VERSION VERSION_LESS 2.8.7 )
	enable_language ( C )
endif ( CMAKE_VERSION VERSION_LESS 2.8.7 ) 

find_package ( OpenMP ${SCAI_FIND_PACKAGE_FLAGS} )

include ( Functions/setAndCheckCache )
setAndCheckCache ( OPENMP )

# LAMA irrelevant entries will be marked as advanced ( Remove them from default cmake GUI )
mark_as_advanced ( OpenMP_C_FLAGS )

## get OpenMP version
if    ( OPENMP_FOUND )
	    try_run ( OPENMP_RUN_RESULT_VAR OPENMP_COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}
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

else  ( OPENMP_FOUND AND USE_OPENMP )

	# Supress unknown pragma warnings if OpenMP is disabled

	if    ( CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES Intel )
        set ( SCAI_WARNING_FLAGS "${SCAI_WARNING_FLAGS} -Wno-unknown-pragmas" )
    endif ( CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES Intel )
    
endif ( OPENMP_FOUND AND USE_OPENMP )

if    ( USE_OPENMP AND NOT OPENMP_FOUND )
    message( FATAL_ERROR "Build of LAMA with OpenMP support enabled, but configuration is incomplete!")
endif ( USE_OPENMP AND NOT OPENMP_FOUND )