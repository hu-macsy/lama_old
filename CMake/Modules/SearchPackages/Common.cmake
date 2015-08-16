###
 # @file SearchPackages.cmake
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
 # @brief List of required and optional packages
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

# Find required packages
set ( REQUIRED_PACKAGES_TO_FIND
        Threads # use ${CMAKE_THREAD_LIBS_INIT} for target_link_libraries
        #add required packages here
    )
    
# Find optional packages
set ( OPTIONAL_PACKAGES_TO_FIND
        #add optional packages here
    )

###  Here we use PThread library for threads
###  Note: FindThreads in CMake is available as Module, but is buggy, needs update of CheckIncludeFiles.cmake
#find_library ( PTHREADS_LIBRARY NAMES pthread pthreads )
