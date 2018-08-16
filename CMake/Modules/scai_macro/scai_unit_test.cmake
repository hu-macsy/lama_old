###
 # @file scai_unit_test.cmake
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
 # @brief Macro to build unit test for SCAI modules
 # @author Thomas Brandes
 # @date 03.07.2017
###

###  Macro to define a unit test executable in a SCAI module
###
###  scai_unit_test( EXECUTABLE <executable> FILES source_file1 source_file2 ... )
###
###  is same as add_executable( <excutable> source_file1 source_file2 .... ) but
###  also adds include directory and links with the unit test library and the module library
###
###  Note: tests will never be installed

macro ( scai_unit_test )

    set ( options CUDA )
    set ( oneValueArgs EXECUTABLE )
    set ( multiValueArgs FILES )

    cmake_parse_arguments ( scai_unit_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( "${scai_unit_test_EXECUTABLE}" STREQUAL "" )
        message ( FATAL_ERROR "scai_unit_test: EXECUTABLE argument not specified" )
    endif ()

    include_directories ( ${BOOST_INCLUDE_DIR} )

    # define test executable, will not be a default built target, but wiht make check

    if ( ${scai_unit_test_CUDA} )
        cuda_add_executable ( ${scai_unit_test_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_unit_test_FILES} )
    else ()
        add_executable ( ${scai_unit_test_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_unit_test_FILES} )
    endif ()

    add_dependencies( tests ${scai_unit_test_EXECUTABLE} )

    if ( WIN32 )
        link_directories ( ${Boost_LIBRARY_DIRS} )
        target_link_libraries ( ${scai_unit_test_EXECUTABLE} ${MODULE_LIBRARY} )
    else ()
        target_link_libraries ( ${scai_unit_test_EXECUTABLE} ${MODULE_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} )
    endif ()

    # run the unit test when 'make test' is called

    add_test( ${scai_unit_test_EXECUTABLE} ${scai_unit_test_EXECUTABLE} )

endmacro ()

macro ( scai_add_test )

    set ( oneValueArgs EXECUTABLE )
    set ( multiValueArgs FILES )

    cmake_parse_arguments ( scai_add_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( "${scai_add_test_EXECUTABLE}" STREQUAL "" )
        message ( FATAL_ERROR "scai_add_test: EXECUTABLE argument not specified" )
    endif ()

    # define test executable, will not be a default built target, but wiht make check

    add_executable ( ${scai_add_test_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_add_test_FILES} )

    add_dependencies( tests ${scai_add_test_EXECUTABLE} )

    target_link_libraries ( ${scai_add_test_EXECUTABLE} ${MODULE_LIBRARY} )

    # test will run by separate script

endmacro ()

