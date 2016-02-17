###
 # @file package/IBVERBS.cmake
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
 # @brief findPackage and configuration of IBVERBS
 # @author Thomas Brandes
 # @date 15.02.2016
 # @since 2.0.0
###

 # - Find ibverbs
 #
 # This module looks for ibverbs support and defines the following values
 #  IBVERBS_FOUND                   TRUE if IBVERBS has been found
 #  IBVERBS_INCLUDE_DIR             the include path for IBVERBS
 #  IBVERBS_LIBRARIES               the library to link against

find_path( IBVERBS_INCLUDE_DIR infiniband/verbs.h
    /usr/local/include
    /usr/include
    $ENV{IBVERBS_INCLUDE_PATH}
)

message( STATUS "IBVERBS_INCLUDE_DIR: ${IBVERBS_INCLUDE_DIR}" )

FIND_LIBRARY( IBVERBS_LIBRARIES ibverbs 
    /usr/local/lib
    /usr/lib
    $ENV{IBVERBS_LIBRARY_PATH}
)

message( STATUS "IBVERBS_LIBRARIES: ${IBVERBS_LIBRARIES}" )

include( FindPackageHandleStandardArgs )

find_package_handle_standard_args( IBVERBS
    DEFAULT_MSG
    IBVERBS_INCLUDE_DIR
    IBVERBS_LIBRARIES
)
