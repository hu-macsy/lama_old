###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2011
 # Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 # for Fraunhofer-Gesellschaft
 #
 # Permission is hereby granted, free of charge, to any person obtaining a copy
 # of this software and associated documentation files (the "Software"), to deal
 # in the Software without restriction, including without limitation the rights
 # to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 # copies of the Software, and to permit persons to whom the Software is
 # furnished to do so, subject to the following conditions:
 #
 # The above copyright notice and this permission notice shall be included in
 # all copies or substantial portions of the Software.
 #
 # THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 # IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 # FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 # AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 # LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 # OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 # SOFTWARE.
 # @endlicense
 #
 # @brief User documentation specification file for CMake
 # @author Kai Buschulte
 # @date 12.09.2012
 # $Id$
###

cmake_minimum_required( VERSION 2.8 )

project ( LAMA_USER_DOC CXX )

if ( NOT LAMA_ROOT AND NOT $ENV{LAMA_ROOT} STREQUAL "" )
    message ( FATAL_ERROR "Please define LAMA_ROOT.")
elseif ( NOT EXISTS ${LAMA_ROOT} )
    message ( FATAL_ERROR "The defined LAMA_ROOT (${LAMA_ROOT}) is not a directory." )
endif ()

if ( EXISTS ${LAMA_ROOT}/include )
    include_directories( ${LAMA_ROOT}/include )
elseif ( EXISTS ${LAMA_ROOT}/src/lama )
    include_directories( ${LAMA_ROOT}/src )
else ()
    message ( WARNING "Found no include subdirectories (include or src) will use include directory LAMA_ROOT ${LAMA_ROOT} " )
    include_directories ( ${LAMA_ROOT} )
endif ()

if ( EXISTS ${LAMA_ROOT}/lib )
    link_directories ( ${LAMA_ROOT}/lib )
elseif ( EXISTS ${LAMA_ROOT}/build )
    link_directories( ${LAMA_ROOT}/build )
else ()
    message ( WARNING "Found no link subdirectories (lib or build) will use link directory LAMA_ROOT ${LAMA_ROOT} " )
    link_directories ( ${LAMA_ROOT} )
endif ()

###
## Compiler Flags
###

### do we need this compiler flags -Wl,-rpath=${LAMA_ROOT}/lib ??? ###

set ( CXX_WARNING_FLAGS "-Wl" )

# CMAKE configuration variable that guarantees adding rpath for installed
# libraries; very useful so that installed library can be used without 
# complex settings of LD_LIBRARY_PATH

set ( CMAKE_SKIP_BUILD_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set ( CMAKE_BUILD_WITH_INSTALL_RPATH FALSE )
set ( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

###
## Add new sources here
###
set ( CXX_SOURCES 
     # tutorial
     tutorial/solver.cpp
     tutorial/simple.cpp
     
     # lecture
     lecture/task0.cpp
     lecture/task1a.cpp
     #lecture/task1b.cpp
     lecture/task2.cpp
     lecture/task2a.cpp
     lecture/task3.cpp
     lecture/task4.cpp
     lecture/task5.cpp
)

# create executeables and link them with lama
foreach ( CXX_SOURCE ${CXX_SOURCES} )
    get_filename_component ( EXE_TARGET ${CXX_SOURCE} NAME_WE )
    message ( STATUS "  Creating executeable ${EXE_TARGET}." )
    
    add_executable ( ${EXE_TARGET} ${CXX_SOURCE} )
    
    # link lama and log4lama
    target_link_libraries ( ${EXE_TARGET} lama log4lama )
    
    # executeables should be created in CMAKE_BINARY_DIR
    set_target_properties ( ${EXE_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
endforeach ()
