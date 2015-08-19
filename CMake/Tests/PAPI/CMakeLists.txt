###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2009-2013
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
 # @brief CMakeLists.txt
 # @author Thomas Brandes
 # @date 25.04.2013
 # @since 1.0.0
###

cmake_minimum_required ( VERSION 2.8 )

project( PapiExample )

# make sure that the default is a RELEASE

if ( NOT CMAKE_BUILD_TYPE )
   set( CMAKE_BUILD_TYPE Release CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE )
endif ( NOT CMAKE_BUILD_TYPE )

# Find PAPI and set flags

set ( PAPI /home/brandes/local/papi-5.0.1 )

find_library( PAPI_LIBRARY
              NAMES papi
              HINTS ${PAPI}/lib )

set ( PAPI_LIBRARIES ${PAPI_LIBRARY} )

include_directories( ${PAPI}/include )

add_executable( example example.cpp )

target_link_libraries ( example ${PAPI_LIBRARIES} )

