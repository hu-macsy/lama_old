###
 # @file CheckC++11.txt
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
 # @brief Check if compiler supports C++11 features.
 # @author Thomas Brandes
 # @date 09.06.2015
 # @since 2.0.0
###

#### Check for -std=c++11

include ( CheckCXXCompilerFlag )

if    ( NOT DEFINED CXX_SUPPORTS_C11 )
    CHECK_CXX_COMPILER_FLAG( -std=c++11 CXX_SUPPORTS_C11 )
endif ( NOT DEFINED CXX_SUPPORTS_C11 )

if    ( CXX_SUPPORTS_C11 )

	if    ( CMAKE_COMPILER_IS_GNUCXX )
        set ( SCAI_LANG_FLAGS "-std=c++11" )
	endif ( CMAKE_COMPILER_IS_GNUCXX )

	if    ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
		set ( SCAI_LANG_FLAGS "-std=c++11" )
	endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

else  ( CXX_SUPPORTS_C11 )
	
	if    ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
		set ( SCAI_LANG_FLAGS "-std=c++0x" )
	endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
	
endif ( CXX_SUPPORTS_C11 )