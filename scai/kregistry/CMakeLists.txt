###
 # @file scai/kregistry/CMakeLists.txt
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
 # @brief scai/kregistry/CMakeLists.txt
 # @author Thomas Brandes
 # @date 30.06.2015
###

## Define name of the SCAI module project

set ( MODULE_NAME kregistry )

set ( INTERNAL_DEPS common logging )
set ( EXTERNAL_DEPS )

## Define all library names with the (global) prefix SCAI_LIBRARY_PREFIX

set ( PROJECT_NAME "${SCAI_LIBRARY_PREFIX}${MODULE_NAME}" )

### find all external packages via the provided wrappers of CMake modules

foreach ( module ${EXTERNAL_DEPS} )
    include( Package/${module} )
endforeach ()

## finalize build flags
include ( ConcludeFlags )

### add definitions

include ( Settings/logging )

add_definitions( ${ADDITIONAL_WARNING_FLAGS} )
add_definitions ( -DSCAI_LOGGING_LEVEL_${SCAI_LOGGING_LEVEL} )
add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

if    ( WIN32 )
    add_definitions ( -DCOMMON_COMPILING_DLL )
endif ( WIN32 )

### Project System-Include ( also configured files )
set( CONFIGURE_INCLUDE_DIR ${CMAKE_BINARY_DIR}/include )
include_directories ( ../.. ${CONFIGURE_INCLUDE_DIR} )

### add includes for external packages

foreach ( module ${EXTERNAL_DEPS} )
    string ( TOUPPER ${module} upper_module )
    include_directories( ${SCAI_${upper_module}_INCLUDE_DIR} )
endforeach ()

### Add classes/headers to source files CXX_SOURCES, CXX_HEADERS, ...

scai_project (

    CLASSES                  # .cpp, .hpp

        ContextFunction
        KernelRegistry

    HEADERS                  # .hpp only

        KernelContextFunction
)

## add subdirectories

if ( BUILD_EXAMPLES )
    add_subdirectory ( examples )
endif ( BUILD_EXAMPLES )

add_subdirectory ( exception )
add_subdirectory ( mepr )

if ( FOUND_BOOST_TEST AND BUILD_TEST )
    add_subdirectory ( test )
endif ( FOUND_BOOST_TEST AND BUILD_TEST )

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

## install library.hpp
install ( FILES ../${MODULE_NAME}.hpp DESTINATION include/scai )

## add custom target doc_${MODULE_NAME} to generate Sphinx user documentation 
include ( CustomCommands/SphinxDoc )
