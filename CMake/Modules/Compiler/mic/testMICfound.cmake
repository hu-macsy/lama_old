###
 # @file testMICfound.cmake
 #
 # @license
 # Copyright (c) 2009-2014
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
 # @brief Detect MIC availability with test program
 # @author Lauretta Schubert
 # @date 04.04.2016
 # @since 2.0.0
###

if ( CMAKE_CXX_COMPILER_ID MATCHES Intel AND USE_OPENMP )

	execute_process ( COMMAND micinfo
					  RESULT_VARIABLE MIC_RUN_RESULT_VAR
					  OUTPUT_VARIABLE MIC_RUN_OUTPUT_VAR)

	# if micinfo return with SUCCESS MIC_RUN_RESULT_VAR is 0
	if    ( NOT MIC_RUN_RESULT_VAR )
		set ( USE_MIC TRUE )
	endif ( NOT MIC_RUN_RESULT_VAR )

else  ( CMAKE_CXX_COMPILER_ID MATCHES Intel AND USE_OPENMP )
	set ( USE_MIC FALSE )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel AND USE_OPENMP )