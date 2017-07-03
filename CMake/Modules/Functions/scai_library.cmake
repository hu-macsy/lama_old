###
 # @file scai_library.cmake
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

##  Be careful: global variables are used
##
##   MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
##   SCAI_LIBRARY_PREFIX, SCAI_VERSION
##   CXX_SOURCES

macro ( scai_library )

    ## Define all library names with the (global) prefix SCAI_LIBRARY_PREFIX

    set ( MODULE_LIBRARY "${SCAI_LIBRARY_PREFIX}${MODULE_NAME}" )

    ### add library ###

    add_library ( ${MODULE_LIBRARY} ${SCAI_LIBRARY_TYPE} ${CXX_SOURCES} )

    set_target_properties ( ${MODULE_LIBRARY} PROPERTIES VERSION ${SCAI_VERSION} )

    ## link internal libraries via internal dependencies, but add library prefix

    string ( REGEX REPLACE "([a-z]+)" "${SCAI_LIBRARY_PREFIX}\\1" INTERNAL_LIBS "${INTERNAL_DEPS}")

    target_link_libraries ( ${MODULE_LIBRARY} ${INTERNAL_LIBS} )

    ## link external libraries via external libraries, SCAI external package wrappers respect naming conventions

    foreach ( module ${EXTERNAL_DEPS} )
        string ( TOUPPER ${module} upper_module )
        target_link_libraries( ${MODULE_LIBRARY} ${SCAI_${upper_module}_LIBRARIES} )
    endforeach ()

    ## install lib

    install ( TARGETS ${MODULE_LIBRARY} DESTINATION lib )

endmacro ()
