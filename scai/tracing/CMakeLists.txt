###
 # @file scai/tracing/CMakeLists.txt
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
 # @brief CMake configuration file for SCAI module tracing
 # @author Thomas Brandes
 # @date 04.07.2017
###

## Define SCAI module project by name and internal/external dependencies

scai_module(

    MODULE_NAME   tracing
    INTERNAL_DEPS common logging
    EXTERNAL_DEPS 
)

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        CallTreeTable
        FileTable
        RegionTable
        TraceConfig
        TraceData
        TraceRegionRecord
        VTInterface

    HEADERS                  # .hpp only

        CallStack
        Counters
        RegionEntry
)

## add subdirectories

scai_subdirectories ( EXAMPLES TEST gui macros )

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ( PREFIX  ${SCAI_LIBRARY_PREFIX}
               TYPE    ${SCAI_LIBRARY_TYPE}
               VERSION ${SCAI_VERSION}
               ${CXX_SOURCES} )

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install library.hpp
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )
