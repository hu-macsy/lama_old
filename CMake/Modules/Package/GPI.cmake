###
 # @file Package/GPI.cmake
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
 # @brief findPackage and configuration of GPI
 # @author Thomas Brandes
 # @date 15.02.2016
###

### USE_GPI              - if GPI is enabled
### GPI_FOUND            - if GPI2 and IBVERBS are founds
### GPI_ENABLED          - if GPI_FOUND AND USE_GPI
### SCAI_GPI_INCLUDE_DIR - GPI include directories (GPI2_INCLUDE_DIR and IBVERBS_INCLUDE_DIR)
### SCAI_GPI_LIBRARIES   - all needed GPI libraries (GPI2_LIBRARIES and IBVERBS_LIBRARIES)
### GPI2_VERSION
### GPI2_INCLUDE_DIR
### IBVERBS_INCLUDE_DIR

set( SCAI_CMAKE_VERBOSE True )

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

set ( GPI_FOUND FALSE )

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

    if    ( SCAI_CMAKE_VERBOSE )
        message( STATUS "GPI2_COMPILE_RESULUT_VAR=${GPI2_COMPILE_RESULT_VAR}" )
        message( STATUS "GPI2_RUN_RESULT_VAR=${GPI2_RUN_RESULT_VAR}" )
        message( STATUS "GPI2_RUN_OUTPUT_VAR=${GPI2_RUN_OUTPUT_VAR}" )
    endif ( SCAI_CMAKE_VERBOSE )

    if ( GPI2_COMPILE_RESULT_VAR )
        # If we could compile correctly we have found GPI 
        # This solution works also if GPI2 has been compiled for Etherent without IBVERBS
        set ( GPI_FOUND TRUE )
    endif ( GPI2_COMPILE_RESULT_VAR )

    set ( GPI2_VERSION ${GPI2_RUN_OUTPUT_VAR} )

endif ( GPI2_FOUND )

### ALLOW to switch off GPI2 explicitly ###
# do what setAndCheckCache does but with 2 packages
# Check if cache variable is already set

if    ( DEFINED USE_GPI )
    # do nothing
    # if cache variable is NOT set
    set( USE_GPI ${USE_GPI} )
else ( DEFINED USE_GPI )
    # Check if package was found
    set ( USE_PACKAGE ${GPI_FOUND} )
    set ( USE_GPI ${USE_PACKAGE} )
endif ( DEFINED USE_GPI )

set ( USE_GPI ${USE_GPI} CACHE BOOL "Enable / Disable use of GPI" )

set ( GPI_ENABLED FALSE )

if    ( USE_GPI AND GPI_FOUND )
    # conclude GPI2 and IBVERBS to SCAI_GPI
    set ( SCAI_GPI_INCLUDE_DIR ${GPI2_INCLUDE_DIR} )
    set ( SCAI_GPI_LIBRARIES ${GPI2_LIBRARIES} ${GPI2_EXTRA_LIBRARIES} )

    set ( GPI_ENABLED TRUE )

endif ( USE_GPI AND GPI_FOUND )
