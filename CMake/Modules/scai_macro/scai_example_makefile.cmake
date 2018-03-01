###
 # @file scai_example_makefile.cmake
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
 # @brief Macros for building examples in SCAI module projects
 # @author Thomas Brandes
 # @date 03.07.2017
###

## This macro generates make.inc and makefile for examples in the installation 
## directory.
##
##  - make.inc is configured from examples_make.inc.in
##  - Makefile is configured from examples_Makefile.in
##
## Note: the configured files are generated in the binary directory before they 
##       will be installed later. Due to name conflics with the cmake build 
##       the file makefile is temporarily renamed to install_Makefile.
##
## The following variables in make.inc/makefile are configured:
##
##   - CMAKE_CXX_COMPILER, CMAKE_CXX_FLAGS, CMAKE_CXX_FLAGS_RELEASE
##   - CMAKE_EXE_LINKER_FLAGS
##   - CUDA_NVCC_EXECUTABLE, CUDA_NVCC_FLAGS_CLEAN
##   - CMAKE_INSTALL_PREFIX
##   - SCAI_BOOST_INCLUDE_DIR
##   - SCAI_CUDA_INCLUDE_DIR  
##   - SCAI_CUDA_LIBRARY_PATH
##   - SCAI_EXAMPLE_LINK_LIBRARIES
##   - SCAI_DEFINES
##   - EXAMPLE_EXECUTABLES
##   - EXAMPLE_MODULES
##   - EXAMPLE_LIBS

macro ( scai_example_makefile )

    set ( options )
    set ( oneValueArgs DESTINATION )
    set ( multiValueArgs )

    cmake_parse_arguments ( scai_example_makefile "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    ## set SCAI_EXAMPLE_LINK_LIBRARIES with own project and all dependent libraries

    set ( SCAI_EXAMPLE_LINK_LIBRARIES "-l${MODULE_LIBRARY}" )

    set ( REVERT_LIST ${INTERNAL_DEPS} ) # because list does not accept variable recursion

    if ( REVERT_LIST ) # is empty for common
        list ( REVERSE REVERT_LIST )
        foreach    ( module ${REVERT_LIST} )
            set ( library "${SCAI_LIBRARY_PREFIX}_${module}" )
            set ( SCAI_EXAMPLE_LINK_LIBRARIES "${SCAI_EXAMPLE_LINK_LIBRARIES} -l${library}" )
        endforeach ( module ${REVERT_LIST} )
    endif ()

    set ( SCAI_EXAMPLE_LINK_LIBRARIES "${SCAI_START_LINK_LIBRARIES} ${SCAI_EXAMPLE_LINK_LIBRARIES} ${SCAI_END_LINK_LIBRARIES}" )

    ## set project specific SCAI_DEFINES

    if ( SCAI_ASSERT_LEVEL )
        set ( SCAI_DEFINES "${SCAI_DEFINES} -DSCAI_ASSERT_LEVEL_${SCAI_ASSERT_LEVEL}" )
    endif ()
    
    if ( SCAI_LOG_LEVEL )
        set ( SCAI_DEFINES "${SCAI_DEFINES} -DSCAI_LOG_LEVEL_${SCAI_LOG_LEVEL}" )
    endif ()

    if ( SCAI_TRACE )
        set ( SCAI_DEFINES "${SCAI_DEFINES} -DSCAI_TRACE_${SCAI_TRACE}" )
    endif ()

    ## SCAI_MODULES = ${EXAMPLE_MODULES} Module.so;DummyModule.so  -> Moudle.so DummyModule.so

    listToString ( " " "${EXAMPLE_EXECUTABLES}" EXAMPLE_EXECUTABLES )
    listToString ( " " "${EXAMPLE_MODULES}" EXAMPLE_MODULES )
    listToString ( " " "${EXAMPLE_LIBS}" EXAMPLE_LIBS )
    listToString ( " " "${CUDA_NVCC_FLAGS}" CUDA_NVCC_FLAGS_CLEAN )
    listToString ( " " "${SCAI_CUDA_LIBRARIES}" SCAI_CUDA_LIBRARIES_CLEAN )

    configure_file ( "${CMAKE_SOURCE_DIR}/examples_make.inc.in" "${CMAKE_CURRENT_BINARY_DIR}/make.inc" @ONLY )
    configure_file ( "${CMAKE_SOURCE_DIR}/examples_Makefile.in" "${CMAKE_CURRENT_BINARY_DIR}/install_Makefile" COPYONLY )

    install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/make.inc 
              DESTINATION ${scai_example_makefile_DESTINATION} )

    install ( FILES       ${CMAKE_CURRENT_BINARY_DIR}/install_Makefile
              DESTINATION ${scai_example_makefile_DESTINATION} 
              RENAME      Makefile  )

endmacro ()

