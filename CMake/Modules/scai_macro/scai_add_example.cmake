###
 # @file scai_add_example.cmake
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
 # @brief Macro to build an exectubale for SCAI modules
 # @author Thomas Brandes
 # @date 03.07.2017
###

###  Macro to define an example executable in a SCAI module
###
###  scai_add_example ( EXECUTABLE <executable> [CUDA] FILES source_file1 source_file2 ... )
###
###  is same as add_executable( <excutable> source_file1 source_file2 .... ) but
###  also adds include directory and links with the module library
###
###  Note: all example executables are only built if target examples is built

macro ( scai_add_example )

    set ( options CUDA )
    set ( oneValueArgs EXECUTABLE LIBRARY )
    set ( multiValueArgs FILES )

    cmake_parse_arguments ( scai_add_example "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    string ( LENGTH "${scai_add_example_EXECUTABLE}" LEN_EXECUTABLE )
    string ( LENGTH "${scai_add_example_LIBRARY}" LEN_LIBRARY )

    # message ( STATUS "add example, LIBRARY ${LEN_LIBRARY} EXECUTABLE ${LEN_EXECUTABLE}" )

    if ( LEN_EXECUTABLE GREATER 0 )

        if ( LEN_LIBRARY GREATER 0 )
            message( FATAL_ERROR "LIBRARY and EXECTUABLE cannot be used together" )
        endif ()
    
        if ( ${scai_add_example_CUDA} )
            cuda_add_executable ( ${scai_add_example_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_add_example_FILES} )
        else ()
            add_executable ( ${scai_add_example_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_add_example_FILES} )
        endif ()

        set ( EXAMPLE_EXECUTABLES ${EXAMPLE_EXECUTABLES} ${scai_add_example_EXECUTABLE} )

        target_link_libraries ( ${scai_add_example_EXECUTABLE} ${MODULE_LIBRARY} )

        # just make sure that all dynamic libraries are linked even if not needed

        set_target_properties ( ${scai_add_example_EXECUTABLE} PROPERTIES LINK_FLAGS "${SCAI_START_LINK_LIBRARIES}"  )

        add_dependencies( examples ${scai_add_example_EXECUTABLE} )

    else ()

        if ( LEN_LIBRARY EQUAL 0 )
            message( FATAL_ERROR "neither LIBRARY nor EXECUTABLE specified" )
        endif ()

        add_library ( ${scai_add_example_LIBRARY} SHARED EXCLUDE_FROM_ALL ${scai_add_example_FILES} )

        set ( EXAMPLE_MODULES ${EXAMPLE_MODULES} "${scai_add_example_LIBRARY}.so" )

        target_link_libraries ( ${scai_add_example_LIBRARY} ${MODULE_LIBRARY} )

        add_dependencies( examples ${scai_add_example_LIBRARY} )

    endif ()

    set ( EXAMPLE_FILES ${EXAMPLE_FILES} ${scai_add_example_FILES} )

endmacro ()

macro ( scai_init_examples )

    set ( EXAMPLE_FILES "" )
    set ( EXAMPLE_EXECUTABLES "" )
    set ( EXAMPLE_MODULES "" ) 

endmacro ()

