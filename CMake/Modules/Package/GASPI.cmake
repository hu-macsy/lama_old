###
 # @file Package/GASPI.cmake
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
 # @brief findPackage and configuration of GASPI
 # @author Thomas Brandes
 # @date 15.02.2016
###

include ( scai_macro/scai_pragma_once )
include ( scai_macro/scai_build_variable )
include ( scai_macro/scai_summary )

### USE_GASPI              - if GASPI is enabled
### GASPI_FOUND            - if GASPI required libraries (e.g. IBVERBS) are found
### SCAI_GASPI_INCLUDE_DIR - GPI include directories (GPI2_INCLUDE_DIR)
### SCAI_GASPI_LIBRARIES   - all needed GPI libraries (GPI2_LIBRARIES and GPI2_EXTRA_LIBRARIES)

# set( SCAI_CMAKE_VERBOSE True )

scai_pragma_once ()

find_package ( GPI2 ${SCAI_FIND_PACKAGE_FLAGS} )

if    ( SCAI_CMAKE_VERBOSE )
    message ( STATUS "GPI2_FOUND=${GPI2_FOUND}" )
    message ( STATUS "GPI2_INCLUDE_DIR=${GPI2_INCLUDE_DIR}" )
    message ( STATUS "GPI2_LIBRARIES=${GPI2_LIBRARIES}" )
endif ( SCAI_CMAKE_VERBOSE )

find_package ( Ibverbs ${SCAI_FIND_PACKAGE_FLAGS} )

if    ( SCAI_CMAKE_VERBOSE )
    message ( STATUS "IBVERBS_FOUND=${IBVERBS_FOUND}" )
    message ( STATUS "IBVERBS_INCLUDE_DIR=${IBVERBS_INCLUDE_DIR}" )
    message ( STATUS "IBVERBS_LIBRARIES=${IBVERBS_LIBRARIES}" )
endif ( SCAI_CMAKE_VERBOSE )

if ( IBVERBS_FOUND )
    set ( GPI2_EXTRA_LIBRARIES ${IBVERBS_LIBRARIES} )
else ( IBVERBS_FOUND )
    set ( GPI2_EXTRA_LIBRARIES "-lpthread" )
endif ( IBVERBS_FOUND )

# ToDo: more convenient analysis which version of GPI-2 (Infiniband/Ethernet, with/without GPU) 
#       has been installed at GPI2_ROOT

set ( GASPI_FOUND FALSE )

if    ( GPI2_FOUND )
    ## get GPI2 version
    try_run ( GPI2_RUN_RESULT_VAR GPI2_COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}/VersionCheck
        ${CMAKE_MODULE_PATH}/VersionCheck/gpi.cpp
        CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${GPI2_INCLUDE_DIR}
        LINK_LIBRARIES ${GPI2_LIBRARIES} ${GPI2_EXTRA_LIBRARIES}
        COMPILE_OUTPUT_VARIABLE GPI2_COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE GPI2_RUN_OUTPUT_VAR )

    if ( SCAI_CMAKE_VERBOSE )
        message( STATUS "GPI2_COMPILE_RESULUT_VAR=${GPI2_COMPILE_RESULT_VAR}" )
        message( STATUS "GPI2_RUN_RESULT_VAR=${GPI2_RUN_RESULT_VAR}" )
        message( STATUS "GPI2_RUN_OUTPUT_VAR=${GPI2_RUN_OUTPUT_VAR}" )
    endif ( SCAI_CMAKE_VERBOSE )

    if ( GPI2_COMPILE_RESULT_VAR )
        # If we could compile correctly we have found GPI 
        # This solution works also if GPI2 has been compiled for Etherent without IBVERBS
        set ( GASPI_FOUND TRUE )
    endif ( GPI2_COMPILE_RESULT_VAR )

    set ( GPI2_VERSION ${GPI2_RUN_OUTPUT_VAR} )

endif ( GPI2_FOUND )

### ALLOW to disable GASPI explicitly ###

scai_build_variable ( NAME      USE_GASPI
                      BOOL 
                      DEFAULT   ${GASPI_FOUND}
                      DOCSTRING "use of GASPI (Global Address Space Programming Interface)" )

if ( USE_GASPI AND GASPI_FOUND )

    # conclude GPI2 and IBVERBS to SCAI_GASPI
    set ( SCAI_GASPI_INCLUDE_DIR ${GPI2_INCLUDE_DIR} )
    set ( SCAI_GASPI_LIBRARIES ${GPI2_LIBRARIES} ${GPI2_EXTRA_LIBRARIES} )

endif ()

scai_summary_external ( NAME      GASPI 
                        ENABLED   ${USE_GASPI}
                        FOUND     ${GASPI_FOUND} 
                        VERSION   "GPI2 ${GPI2_VERSION}"
                        INCLUDE   ${SCAI_GASPI_INCLUDE_DIR} 
                        LIBRARIES ${SCAI_GASPI_LIBRARIES} )

