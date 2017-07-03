###
 # @file scai/dmemo/CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2017
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
 # @brief CMake specification file for project dmemo (distributed memory)
 # @author Thomas Brandes
 # @date 16.02.2016
###

## Define name of the SCAI module project

scai_module(

    MODULE_NAME   dmemo
    INTERNAL_DEPS common logging tracing tasking hmemo kregistry blaskernel utilskernel sparsekernel 
    EXTERNAL_DEPS OpenMP MPI GPI GraphPartitioning 
)

## MetisDistribution only added as source if graph partitioning is available

if    ( ( METIS_FOUND OR PARMETIS_FOUND ) AND USE_GRAPHPARTITIONING )
    if    ( SCAI_INDEX_TYPE STREQUAL "int" )
        set ( METISDISTRIBUTION MetisDistribution )
    endif ( SCAI_INDEX_TYPE STREQUAL "int" )
endif ( ( METIS_FOUND OR PARMETIS_FOUND ) AND USE_GRAPHPARTITIONING )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        Communicator
        NoCommunicator
        Distributed
        CommunicationPlan

        BlockDistribution
        CyclicDistribution
        Distribution
        GenBlockDistribution
        GeneralDistribution
        GridDistribution
        Halo
        HaloBuilder
        NoDistribution
        Redistributor
        SingleDistribution

        ${METISDISTRIBUTION}

    )

### add subdirectories

if ( GPI_FOUND AND USE_GPI )
    include( Package/Boost )
    message( STATUS "add BOOST_INCLUDE_DIR=${BOOST_INCLUDE_DIR}" )
    include_directories ( ${BOOST_INCLUDE_DIR} )
    add_subdirectory( gpi )
endif ( GPI_FOUND AND USE_GPI )

if ( MPI_FOUND AND USE_MPI )
    add_subdirectory ( mpi )
endif ( MPI_FOUND AND USE_MPI )

scai_subdirectories( EXAMPLES TEST )

### set link directories ###
if ( WIN32 )
    link_directories ( ${Boost_LIBRARY_DIRS} )
endif ( WIN32 )

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ()

### install ###

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install "<library>.hpp"
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )

