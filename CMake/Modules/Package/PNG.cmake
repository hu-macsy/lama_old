###
 # @file Package/PNG.cmake
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
 # @brief SCAI wrapper for find_package( PNG )
 # @author Thomas Brandes
 # @date 14.05.2017
###

include ( scai_macro/scai_pragma_once )

scai_pragma_once()

#  CMake provides already a module to find the PNG reference library, use it

find_package( PNG ${SCAI_FIND_PACKAGE_FLAGS} )

# returns PNG_FOUND, PNG_INCLUDE_DIRS (cache), PNG_LIBRARIES

# message( STATUS "PNG_FOUND=${PNG_FOUND}" )
# message( STATUS "PNG_INCLUDE_DIRS=${PNG_INCLUDE_DIRS}" )
# message( STATUS "PNG_LIBRARIES=${PNG_LIBRARIES}" )

# Now make some adaptions to fit these variables to the SCAI project rules

scai_build_variable ( NAME      USE_PNG
                      BOOL 
                      DEFAULT   ${PNG_FOUND}
                      DOCSTRING "use of PNG library (read/write PNG images)" )

# set the corresponding SCAI variables to inherit automatic settings by external dependencies

if ( PNG_FOUND )

    #  set SCAI_PNG_xxx variable, in CACHE so they might be used by any module

    set ( SCAI_PNG_LIBRARIES ${PNG_LIBRARIES} CACHE PATH "PNG library" )
    set ( SCAI_PNG_INCLUDE_DIR ${PNG_INCLUDE_DIRS} CACHE PATH "PNG include directory" )

    ## get PNG version

    try_run ( PNG_RUN_RESULT_VAR PNG_COMPILE_RESULT_VAR
        ${CMAKE_BINARY_DIR}/VersionCheck
        ${CMAKE_MODULE_PATH}/VersionCheck/pnglib.cpp
        CMAKE_FLAGS 
        -DINCLUDE_DIRECTORIES:STRING=${PNG_INCLUDE_DIRS}
        COMPILE_OUTPUT_VARIABLE PNG_COMPILE_OUTPUT_VAR
        RUN_OUTPUT_VARIABLE PNG_RUN_OUTPUT_VAR )

    set ( PNG_VERSION ${PNG_RUN_OUTPUT_VAR} )

elseif ( SCAI_USE_PNG )

    message ( ERROR "PNG not found" )

endif ()

mark_as_advanced ( SCAI_PNG_LIBRARIES )
mark_as_advanced ( SCAI_PNG_INCLUDE_DIR )

scai_summary_external ( NAME      "PNG (image library)"
                        ENABLED   ${USE_PNG}
                        FOUND     ${PNG_FOUND} 
                        VERSION   ${PNG_VERSION} 
                        INCLUDE   ${PNG_INCLUDE_DIR} 
                        LIBRARIES ${PNG_LIBRARIES}    )

set ( SCAI_PNG_CHECK_DONE True )
