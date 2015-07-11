###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2015
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
 # @brief CMake configuration file for logging examples
 # @author Thomas Brandes
 # @date 16.06.2015
###

# Try to find OpenMP compile and link flags for the OpenMP example

# LogTest example, here we compile also LOG_TRACE macros

add_executable( LogTest LogTest.cpp )
target_link_libraries( LogTest ${logging_LIBRARY} )
set_target_properties( LogTest PROPERTIES COMPILE_FLAGS -DLAMA_LOG_LEVEL_TRACE )

# ThreadTest example, here we need the Boost Thread library

# Boost thread required

set ( Boost_DEBUG TRUE )
set ( Boost_DETAILED_FAILURE_MSG TRUE )

find_package ( Boost COMPONENTS thread system )

if ( Boost_THREAD_FOUND )
    if ( Boost_SYSTEM_FOUND )
        add_executable( ThreadTest ThreadTest.cpp )
        target_link_libraries( ThreadTest ${logging_LIBRARY} ${Boost_THREAD_LIBRARY} ${Boost_SYSTEM_LIBRARY} )
    else ()
        message ( STATUS "Boost system library not found" )
    endif ()
else ()
    message ( STATUS "Boost thread library not found" )
endif ()

# OpenMPTest example

find_package ( OpenMP )

if ( OPENMP_FOUND )
    # Adding this example only possible if OpenMP is available
    add_executable( OpenMPTest OpenMPTest.cpp )
    target_link_libraries( OpenMPTest ${logging_LIBRARY} )
    set_target_properties( OpenMPTest PROPERTIES COMPILE_FLAGS ${OpenMP_CXX_FLAGS} )
    set_target_properties( OpenMPTest PROPERTIES LINK_FLAGS ${OpenMP_CXX_FLAGS} )
else ()
    message ( STATUS "OpenMP example not built, no compiler option found" )
endif()

