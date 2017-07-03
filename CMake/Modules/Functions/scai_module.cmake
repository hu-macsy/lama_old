###
 # @file scai_module.cmake
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
 # @brief CMake macro to define a SCAI module library                          
 # @author Thomas Brandes
 # @date 04.07.2017
###

##  Macro to define a SCAI module project
##
##  Output variables (for further use):
##
##    MODULE_NAME
##    PROJECT_NAME
##    INTERNAL_DEPS
##    EXTERNAL_DEPS
###

macro ( scai_module )

    # specify the keywors supported in arguments 

    set ( options )
    set ( oneValueArgs MODULE_NAME )
    set ( multiValueArgs INTERNAL_DEPS EXTERNAL_DEPS )

    cmake_parse_arguments ( scai_module "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    # message ( STATUS "scai_module MODULE_NAME = ${scai_module_MODULE_NAME}" )
    # message ( STATUS "scai_module INTERNAL_DEPS = ${scai_module_INTERNAL_DEPS}" )
    # message ( STATUS "scai_module EXTERNAL_DEPS = ${scai_module_EXTERNAL_DEPS}" )

    ## Define all library names with the (global) prefix SCAI_LIBRARY_PREFIX

    set ( MODULE_NAME "${scai_module_MODULE_NAME}" )
    set ( PROJECT_NAME "${SCAI_LIBRARY_PREFIX}${MODULE_NAME}" )
  
    set ( INTERNAL_DEPS ${scai_module_INTERNAL_DEPS} )
    set ( EXTERNAL_DEPS ${scai_module_EXTERNAL_DEPS} )

    ### find all external packages via the provided wrappers of CMake modules

    foreach ( module ${EXTERNAL_DEPS} )
        include( Package/${module} )
    endforeach ()

    add_definitions ( ${ADDITIONAL_WARNING_FLAGS} )
    add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

    foreach ( module ${INTERNAL_DEPS} )

        if ( ${module} STREQUAL "logging" )
            include ( Settings/logging )
            add_definitions ( -DSCAI_LOGGING_LEVEL_${SCAI_LOGGING_LEVEL} )
        endif ()

        if ( ${module} STREQUAL "tracing" )
            include ( Settings/tracing )
            add_definitions ( -DSCAI_TRACING_${SCAI_TRACING} )
        endif ()

    endforeach ( module ${INTERNAL_DEPS} )
    
    if ( WIN32 )
        add_definitions ( -DCOMMON_COMPILING_DLL )
    endif ( WIN32 )

    ### set include directories ( same for all module projects ) 

    include_directories ( ${CMAKE_SOURCE_DIR}/.. )       #  for all internal projects
    include_directories ( ${CMAKE_BINARY_DIR}/include )  #  for all configured includes

    ### add includes for external packages

    foreach ( module ${EXTERNAL_DEPS} )
        string ( TOUPPER ${module} upper_module )
        include_directories( ${SCAI_${upper_module}_INCLUDE_DIR} )
    endforeach ()

endmacro ()
