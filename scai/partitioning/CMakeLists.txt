###
 # @file scai/partitioning/CMakeLists.txt
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
 # @brief CMake configuration file for SCAI module partitioning
 # @author Thomas Brandes
 # @date 23.08.2017
###

## Include required scai macros

include ( scai_macro/scai_module )
include ( scai_macro/scai_project )
include ( scai_macro/scai_subdirectories )
include ( scai_macro/scai_library )
include ( scai_function/relative_install )

## Define new SCAI module project with its internal/external dependencies

scai_module(

    MODULE_NAME   partitioning
    INTERNAL_DEPS common logging tracing hmemo kregistry utilskernel dmemo lama
    EXTERNAL_DEPS ParMetis MPI Metis
)

# ToDo: The order ParMetis < Metis is important here as ParMetis links also with Metis

set ( PARTITIONING_CLASSES Partitioning BlockPartitioning CyclicPartitioning CSRGraph CSRGraph2 )

if ( USE_METIS )
    set ( PARTITIONING_CLASSES ${PARTITIONING_CLASSES} MetisPartitioning )
endif ()
    
# ParMetis depends on Metis and MPI

if ( USE_PARMETIS AND USE_METIS AND USE_MPI )
    set ( PARTITIONING_CLASSES ${PARTITIONING_CLASSES} ParMetisPartitioning )
endif ()

# message ( STATUS "PARTITIONING_CLASSES=${PARTITIONING_CLASSES}" )

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        ${PARTITIONING_CLASSES}
    )

scai_subdirectories( EXAMPLES TEST )

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ( PREFIX  ${SCAI_LIBRARY_PREFIX}
               TYPE    ${SCAI_LIBRARY_TYPE}
               VERSION ${SCAI_VERSION}
               ${CXX_SOURCES} )

### install ###

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install "<library>.hpp"
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )

