###
 # @file scai_test.cmake
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
 # @brief Macro to test programs for SCAI modules
 # @author Thomas Brandes
 # @date 03.07.2017
###

###  Macro to define a test executable in a SCAI module
###
###  scai_test( EXECUTABLE <executable> FILES source_file1 source_file2 ... )
###
###   - same as add_executable( <excutable> source_file1 source_file2 .... ) but
###   - executable will be built only if targets tests (or check) is specified
###   - always links with the 'current' module library
###
###  Further options:
###
###   - CUDA      if set executable is built via cuda_add_executable
###   - RUN       if set executable runs with make test
###   - UNIT_TEST test is built with boost unit test framework
###
###  Note: tests will never be installed

macro ( scai_test )

    set ( options CUDA UNIT_TEST RUN )
    set ( oneValueArgs EXECUTABLE )
    set ( multiValueArgs FILES )

    cmake_parse_arguments ( scai_test "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    if ( "${scai_test_EXECUTABLE}" STREQUAL "" )
        message ( FATAL_ERROR "scai_test: EXECUTABLE argument not specified" )
    endif ()

    # define test executable, will not be a default built target, but wiht make check

    if ( ${scai_test_CUDA} )
        cuda_add_executable ( ${scai_test_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_test_FILES} )
    else ()
        add_executable ( ${scai_test_EXECUTABLE} EXCLUDE_FROM_ALL ${scai_test_FILES} )
    endif ()

    add_dependencies( tests ${scai_test_EXECUTABLE} )

    target_link_libraries ( ${scai_test_EXECUTABLE} ${MODULE_LIBRARY} )

    if ( ${scai_test_UNIT_TEST} )
        if (NOT TESTSUPPORT_INCLUDE_DIR)
            message(FATAL_ERROR "TESTSUPPORT_INCLUDE_DIR is not set, can not correctly prepare test executable for compilation.")
        endif()

        target_include_directories( ${scai_test_EXECUTABLE} SYSTEM PRIVATE ${Boost_INCLUDE_DIRS} ${TESTSUPPORT_INCLUDE_DIR})

        if ( WIN32 )
            link_directories ( ${Boost_LIBRARY_DIRS} )
        else ()
            target_link_libraries ( ${scai_test_EXECUTABLE} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} testsupport)
        endif ()
    endif ()

    if ( ${scai_test_RUN} )

        # run the unit test when 'make test' is called

        add_test( ${scai_test_EXECUTABLE} ${scai_test_EXECUTABLE} )

    endif()

endmacro ()

