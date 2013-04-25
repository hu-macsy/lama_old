 ###
 # @file FindGASNET.cmake
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
 # @brief Find GASNET
 # @author
 # @date 25.04.2013
###

 # - Find GASNET
 #
 # This module looks for GASNET support and defines the following values
 #  GASNET_FOUND                   TRUE if GASNET has been found
 #  GASNET_LIBRARIES               the library to link against
   
#set( GASNET_INCLUDE_DIR "/usr/local/lib" ) 
   
FIND_LIBRARY( GASNET_LIBRARIES gasnet-ibv-par
	/usr/local/lib
	/usr/lib
	$ENV{GASNET_LIBRARY_PATH}
)
   
INCLUDE( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( GASNET
    DEFAULT_MSG
    GASNET_LIBRARIES
)

MARK_AS_ADVANCED( GASNET_LIBRARIES )

