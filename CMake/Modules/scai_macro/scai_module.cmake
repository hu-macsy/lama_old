###
 # @file scai_module.cmake
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
 # @brief CMake macro to define a SCAI module with its dependencies
 # @author Thomas Brandes
 # @date 04.07.2017
###

##  Macro to define a SCAI module project
##
##   scai_module ( MODULE_NAME   mymodule 
##                 INTERNAL_DEPS common logging hmemo
##                 EXTERNAL_DEPS BLAS CUDA OpenMP      
##               )
##
##    @param mymodule name of the module (must be unique)
##    @param INTERNAL_DEPS <module_1> ... are names of already defined modules
##    @param EXTERNAL_DEPS <package_1> ... <package_n> are names of external modules
##
##    - for each external PPP a CMake file package/PPP must be available/##
##  Keyword variables remain set for further use in the current scope
##
##    MODULE_NAME
##    INTERNAL_DEPS
##    EXTERNAL_DEPS
##
##   This macro has the following purposes
##
##    - each module of INTERNAL_DEPS is checked whether it has been defined
##    - each package of EXTERNAL_DEPS is configured
##    - definitions of INTERNAL_DEPS and EXTERNAL_DEPS are 
##    - required include files of EXTERNAL_DEPS are set via include_directories
## 

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

    set ( MODULE_NAME    "${scai_module_MODULE_NAME}" )
    set ( MODULE_LIBRARY "${SCAI_LIBRARY_PREFIX}_${MODULE_NAME}" )

    message ( STATUS "Configuring module ${MODULE_NAME} -> builds lib ${MODULE_LIBRARY}" )
  
    set ( INTERNAL_DEPS ${scai_module_INTERNAL_DEPS} )
    set ( EXTERNAL_DEPS ${scai_module_EXTERNAL_DEPS} )

    ### verify that all internal modules have been defined before

    foreach ( module ${INTERNAL_DEPS} )

        list ( FIND SCAI_DEFINED_MODULES ${module}  _index)

        if ( ${_index} EQUAL -1 )

            list ( FIND SCAI_ALL_MODULES ${module} _index )
 
            if ( ${_index} EQUAL -1 )

               set ( EXPLAIN "is not a SCAI module at all" )

            else () 

               set ( EXPLAIN "must be defined before" )

            endif ()

            # message ( FATAL_ERROR "scai_module ${MODULE_NAME}: ${module} of INTERNAL_DEPS not defined, ${EXPLAIN}" )
 
            add_subdirectory ( ../${module} ${CMAKE_BINARY_DIR}/${module} )

        endif ( )
    endforeach ()

    ### find all external packages via the provided wrappers of CMake modules

    foreach ( module ${EXTERNAL_DEPS} )
        include( Package/${module} )
    endforeach ()

    add_definitions ( ${ADDITIONAL_WARNING_FLAGS} )
    add_definitions ( -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL} )

    foreach ( module ${INTERNAL_DEPS} )

        if ( ${module} STREQUAL "logging" )
            add_definitions ( -DSCAI_LOG_LEVEL_${SCAI_LOG_LEVEL} )
        endif ()

        if ( ${module} STREQUAL "tracing" )
            add_definitions ( -DSCAI_TRACE_${SCAI_TRACE} )
        endif ()

    endforeach ( module ${INTERNAL_DEPS} )
    
    if ( WIN32 )
        add_definitions ( -DCOMMON_COMPILING_DLL )
    endif ( WIN32 )

    ### set one include directoriy for all internal SCAI projects

    include_directories ( ${CMAKE_SOURCE_DIR}/.. )       

    ##  set one include directory for all configured includes
    ##  Note: define it as SYSTEM to avoid reconfiguration with each call of make

    include_directories ( SYSTEM ${CMAKE_BINARY_DIR}/include )  #  for all configured includes

    ### add includes for external packages

    foreach ( module ${EXTERNAL_DEPS} )
        string ( TOUPPER ${module} upper_module )
        include_directories( ${SCAI_${upper_module}_INCLUDE_DIR} )
    endforeach ()

    ### add this new module to 'global' list of defined modules

    list( APPEND SCAI_DEFINED_MODULES ${MODULE_NAME} )
    set ( SCAI_DEFINED_MODULES ${SCAI_DEFINED_MODULES} PARENT_SCOPE )

endmacro ()
