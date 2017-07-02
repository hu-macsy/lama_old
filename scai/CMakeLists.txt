###
 # @file CMakeLists.txt.tb
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
 # @brief CMake configuration file for involved SCAI projects.
 # @author Thomas Brandes
 # @date 01.07.2017
###

cmake_minimum_required ( VERSION 2.8.8 )

project ( LAMA_ALL CXX )

### Define all modules that are used ( inclusive all dependent ones )

set ( SCAI_PROJECT_LIST common logging tracing tasking kregistry hmemo blaskernel utilskernel sparsekernel dmemo lama solver )

### Define some general variables used for identifaction 

set ( SCAI_VERSION_NAME "Barrancas Blancas" )
set ( SCAI_COPYRIGHT_YEAR 2017 )
set ( SCAI_VERSION 2.2.0 )
set ( SCAI_LIBRARY_PREFIX "scai_" )             # more unique names for all SCAI module libraries

### Set Module Path

set ( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/../CMake/Modules" )

### set cmake quiet for untyped messages
if    ( NOT SCAI_CMAKE_VERBOSE )
    include ( Functions/setCMakeQuiet )
endif ( NOT SCAI_CMAKE_VERBOSE )

### specify the list of SCAI projects that are used
### dependencies are not resolved here but will be verified during configuration


foreach ( SCAI_PROJECT ${SCAI_PROJECT_LIST} )

    message( STATUS "Configuring SCAI project ${SCAI_PROJECT}" )

    # Note: names of subdirectories are same as SCAI_PROJECT names

    add_subdirectory( ${SCAI_PROJECT} )

endforeach ()

include ( Package/Sphinx )
include ( Package/Doxygen )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 

set ( MODULE_NAME "lama_all" )
set ( INTERNAL_DEPS ${SCAI_PROJECT_LIST} )

include ( CustomCommands/SphinxDoc )

## add custom target doxygendoc

include ( CustomCommands/DoxygenDoc )

## add custom target distClean

include ( CustomCommands/DistClean )

## summary of what is found

## ToDo: Summary uses info from module projects no more available here, so redo it

set ( EXTERNAL_DEPS CUDA MIC dl Thread Java SCAI_BLAS MPI GPI GraphPartitioning ZLIB PNG )

foreach ( package ${EXTERNAL_DEPS} )
    include ( Package/${package} )
endforeach ()

include ( Summaries/scai )

