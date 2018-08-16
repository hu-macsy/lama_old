###
 # @file scai_test_scripts.cmake
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
 # @brief Macros for building test and run scripts
 # @author Thomas Brandes
 # @date 03.07.2017
###

##   scai_test_scripts( SCRITPS       <script_file1> ...
##                      FILES         <file1> <file2> ...
##                      [CONFIGURE]                
##
##   - scripts and files are copied from current source directory to current build directory
##   - scripts get an executable permission
##   - if CONFIGURE is set scripts and files are configured from file with suffix ".in"
##   - if CODE_COVERAGE is true, code_coverage.sh will be configured as script

set ( EXECUTABLE_FILE_PERMISSIONS WORLD_READ WORLD_EXECUTE OWNER_READ GROUP_READ GROUP_EXECUTE OWNER_WRITE OWNER_EXECUTE )

macro ( scai_test_scripts )

    set ( options CONFIGURE )
    set ( multiValueArgs FILES SCRIPTS )

    cmake_parse_arguments ( scai_test_scripts "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )

    ### Copy the scripts in build directory ### 

    foreach ( test_file ${scai_test_scripts_SCRIPTS} )

        if ( ${scai_test_scripts_CONFIGURE} )

            configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/${test_file}.in" 
                            "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${test_file}"
                             @ONLY )

            file ( COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/${test_file} 
                   DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
                   FILE_PERMISSIONS ${EXECUTABLE_FILE_PERMISSIONS} )

        else()

            file ( COPY "${CMAKE_CURRENT_SOURCE_DIR}/${test_file}" 
                   DESTINATION ${CMAKE_CURRENT_BINARY_DIR} 
                   FILE_PERMISSIONS ${EXECUTABLE_FILE_PERMISSIONS} 
                 ) 

        endif()

    endforeach() 

    ### Copy the files in build directory, no permissions set  ### 

    foreach ( test_file ${scai_test_scripts_FILES} )

        if ( ${scai_test_scripts_CONFIGURE} )

            ## configure_file is always from CURRENT_SOURCE to CURRENT_BINARY 

            configure_file ( "${test_file}.in" "${test_file}" )

        else()

            file ( COPY "${CMAKE_CURRENT_SOURCE_DIR}/${test_file}" 
                   DESTINATION ${CMAKE_CURRENT_BINARY_DIR}  )

        endif()

    endforeach() 

endmacro ()

##
##  scai_test_code_coverage( ${USE_CODE_COVERAGE} )
##
##   - generates/copies necessary scripts for code coverage
##

macro ( scai_test_code_coverage FLAG )

    if ( ${FLAG} )

        configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/code_coverage.sh.in" "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/code_coverage.sh" @ONLY)

        file ( COPY ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/code_coverage.sh 
               DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
               FILE_PERMISSIONS ${EXECUTABLE_FILE_PERMISSIONS} )

        file ( COPY ${CMAKE_SOURCE_DIR}/scai_code_coverage_functions.sh 
               DESTINATION ${CMAKE_CURRENT_BINARY_DIR}
               FILE_PERMISSIONS ${EXECUTABLE_FILE_PERMISSIONS} )

    endif ()

endmacro ()

