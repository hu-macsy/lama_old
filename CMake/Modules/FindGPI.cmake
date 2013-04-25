 ###
 # @file FindGPI.cmake
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
 # @brief Find GPI
 # @author
 # @date 25.04.2013
###

 # - Find GPI
 #
 # This module looks for GPI support and defines the following values
 #  GPI_FOUND                   TRUE if GPI has been found
 #  GPI_INCLUDE_DIR             the include path for GPI
 #  GPI_LIBRARIES               the library to link against

FIND_PATH(GPI_INCLUDE_DIR GPI.h
	/usr/local/include
	/usr/include
	$ENV{GPI_INCLUDE_PATH}
)
message( STATUS "GPI_INCLUDE_DIR: ${GPI_INCLUDE_DIR}" )
   
FIND_LIBRARY( GPI_LIBRARIES GPI 
	/usr/local/lib
	/usr/lib
	$ENV{GPI_LIBRARY_PATH}
)

message( STATUS "GPI_LIBRARIES: ${GPI_LIBRARIES}" )
   
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GPI
    DEFAULT_MSG
    GPI_INCLUDE_DIR
    GPI_LIBRARIES
)

MARK_AS_ADVANCED( GPI_INCLUDE_DIR GPI_LIBRARIES )
