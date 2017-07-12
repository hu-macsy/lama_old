###
 # @file Package/ZLIB.cmake
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
 # @brief find ZLIB include path and library
 # @author Thomas Brandes
 # @date 25.11.2016
###

find_package( ZLIB ${SCAI_FIND_PACKAGE_FLAGS} )

# returns ZLIB_FOUND, ZLIB_INCLUDE_DIR (cache), ZLIB_LIBRARY

scai_build_variable ( NAME      USE_ZLIB
                      BOOL 
                      DEFAULT   ${ZLIB_FOUND}
                      DOCSTRING "use of ZLIB libray (data compression)" )

# set the corresponding SCAI variables to inherit automatic settings by external dependencies

if ( ZLIB_FOUND )
    set ( SCAI_ZLIB_LIBRARIES ${ZLIB_LIBRARY} CACHE PATH "ZLIB library" )
    set ( SCAI_ZLIB_INCLUDE_DIR ${ZLIB_INCLUDE_DIR} CACHE PATH "ZLIB include directory" )

    ## get ZLIB version
    try_run ( ZLIB_RUN_RESULT_VAR ZLIB_COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}/VersionCheck
        ${CMAKE_MODULE_PATH}/VersionCheck/zlib.cpp
        CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${ZLIB_INCLUDE_DIR}
        COMPILE_OUTPUT_VARIABLE ZLIB_COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE ZLIB_RUN_OUTPUT_VAR )

    set ( ZLIB_VERSION ${ZLIB_RUN_OUTPUT_VAR} )

elseif( USE_ZLIB )

    message ( ERROR "ZLIB not found" )

endif ()

mark_as_advanced ( SCAI_ZLIB_LIBRARIES )
mark_as_advanced ( SCAI_ZLIB_INCLUDE_DIR )

scai_summary_enabled ( NAME "ZLIB (data compression)" ENABLED ${USE_ZLIB} )

scai_summary_external ( NAME      zlib
                        FOUND     ${ZLIB_FOUND} 
                        VERSION   ${ZLIB_VERSION} 
                        INCLUDE   ${ZLIB_INCLUDE_DIR} 
                        LIBRARIES ${ZLIB_LIBRARY}  )

