###
 # @file scai/tasking/CMakeLists.txt
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
 # @brief CMake configuration file for SCAI module tasking
 # @author Thomas Brandes
 # @date 04.07.2017
###

## Include required scai macros

include ( scai_macro/scai_module )
include ( scai_macro/scai_project )
include ( scai_macro/scai_subdirectories )
include ( scai_macro/scai_library )
include ( scai_function/relative_install )

## Define new SCAI module project with its internal/external dependencies

scai_module (

    MODULE_NAME   tasking
    INTERNAL_DEPS common logging tracing
    EXTERNAL_DEPS CUDA
)

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    HEADERS                  # .hpp only

        NoSyncToken

    CLASSES                  # .cpp, .hpp

        SyncToken
        TaskSyncToken
        Task
        ThreadPool
)

## add subdirectories

scai_subdirectories( CUDA EXAMPLES TEST )

if ( CUDA_FOUND AND USE_CUDA )
    cuda_compile ( CUDA_FILES ${CUDA_SOURCES} )
    set ( CXX_SOURCES ${CXX_SOURCES} ${CUDA_FILES} )
endif()

## Define library via macro using MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
## and global variables SCAI_LIBRARY_PREFIX, SCAI_VERSION, ....

scai_library ( PREFIX  ${SCAI_LIBRARY_PREFIX}
               TYPE    ${SCAI_LIBRARY_TYPE}
               VERSION ${SCAI_VERSION}
               ${CXX_SOURCES} )

### install ###

## install headers
relative_install ( FILES ${CXX_HEADERS} DESTINATION "include/scai/${MODULE_NAME}" )

## install library.hpp
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )

