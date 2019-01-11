###
 # @file scai_library.cmake
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
 # @brief CMake macro to define a SCAI module library                          
 # @author Thomas Brandes
 # @date 04.07.2017
###

##  Be careful: global variables are used
##
##   MODULE_NAME, INTERNAL_DEPS, EXTERNAL_DEPS
##
##   scai_libary ( PREFIX scai_ VERSION 1.0.1 TYPE STATIC <source_files> )
##
##   sets: MODULE_LIBRARY ( = <PREFIX><MODULE_NAME> )

macro ( scai_library )

    # specify the keywors supported in arguments 

    set ( options )
    set ( oneValueArgs PREFIX VERSION TYPE )
    set ( multiValueArgs CLASSES HEADERS SOURCES CUDA_CLASSES CUDA_SOURCES )

    cmake_parse_arguments ( scai_library "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( ${scai_library_TYPE} STREQUAL "SHARED" )
    elseif ( ${scai_library_TYPE} STREQUAL "STATIC" )
    else ()
        message ( ERROR "library TYPE ${scai_library_TYPE} illegal, must be SHARED or STATIC" )
    endif ()

    ## Use for all library names the specified prefix

    set ( MODULE_LIBRARY "${scai_library_PREFIX}_${MODULE_NAME}" )

    ### add library ###

    add_library ( ${MODULE_LIBRARY} ${scai_library_TYPE} ${scai_library_UNPARSED_ARGUMENTS} )

    set_target_properties ( ${MODULE_LIBRARY} PROPERTIES VERSION ${scai_library_VERSION} )

    ## link internal libraries via internal dependencies, but add library prefix

    string ( REGEX REPLACE "([a-z]+)" "${scai_library_PREFIX}_\\1" INTERNAL_LIBS "${INTERNAL_DEPS}")

    list ( REVERSE INTERNAL_LIBS )

    ## LAMA relies on registrations of kernels or registration in a factory during static initialization
    ## so some compilers require  -Wl,--no-as-needed <internal libs> -Wl,--as-needed 

    target_link_libraries ( ${MODULE_LIBRARY} ${SCAI_START_LINK_LIBRARIES} ${INTERNAL_LIBS} ${SCAI_END_LINK_LIBRARIES} )

    ## link external libraries via external libraries, SCAI external package wrappers respect naming conventions

    foreach ( module ${EXTERNAL_DEPS} )
        string ( TOUPPER ${module} upper_module )
        target_link_libraries( ${MODULE_LIBRARY} ${SCAI_${upper_module}_LIBRARIES} )
    endforeach ()

    ## install lib

    install ( TARGETS ${MODULE_LIBRARY} DESTINATION lib )

endmacro ()
