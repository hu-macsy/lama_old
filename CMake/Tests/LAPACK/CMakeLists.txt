###
 # @file CMakeLists.txt
 #
 # @license
 # Copyright (c) 2013
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
 # $Id$
###

cmake_minimum_required (VERSION 2.8 )

# Attention: project( GramSchmidt CXX ) does not work

project( GramSchmidt )

# make sure that the default is a RELEASE

if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

set ( CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../Modules" )

# Example of setting BLAS and LAPACK library for Windows
# 
# set( BLAS_blas_LIBRARY "E:/lapack/libblas.lib" )
# set( LAPACK_lapack_LIBRARY "E:/lapack/liblapack.lib" )

include( FindLAPACK )

if ( LAPACK_FOUND )
    message( STATUS "LAPACK found, LAPACK_LIBRARIES=${LAPACK_LIBRARIES}" )
    add_executable( QR QR.cpp Walltime.cpp )
    target_link_libraries( QR ${LAPACK_LIBRARIES} )
else( LAPACK_FOUND )
    message( STATUS "LAPACK not found" )
endif( LAPACK_FOUND )

