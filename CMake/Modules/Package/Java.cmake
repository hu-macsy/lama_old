###
 # @file Java.cmake
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
 # @brief findPackage and configuration of Java
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

if    ( NOT APPLE )
    if    ( CMAKE_VERSION VERSION_GREATER 2.8.11 )
	# Enable Java Compilation
	# this will set at least CMAKE_Java_COMPILER CMAKE_Java_ARCHIVE 
        find_package( Java )

        # LAMA irrelevant entries will be removed from cmake GUI completely
        set ( Java_JAR_EXECUTABLE "${Java_JAR_EXECUTABLE}" CACHE INTERNAL "" )
        set ( Java_JAVADOC_EXECUTABLE "${Java_JAVADOC_EXECUTABLE}" CACHE INTERNAL "" )
        set ( Java_JAVAH_EXECUTABLE "${Java_JAVAH_EXECUTABLE}" CACHE INTERNAL "" )
        set ( Java_JAVA_EXECUTABLE "${Java_JAVA_EXECUTABLE}" CACHE INTERNAL "" )
    endif ( CMAKE_VERSION VERSION_GREATER 2.8.11 )
endif ( NOT APPLE )
