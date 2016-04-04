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


if ( CMAKE_CXX_COMPILER_ID MATCHES Intel )

    # MIC_COMPILE_RESULT_VAR is TRUE when compile succeeds
    # MIC_RUN_RESULT_VAR is zero when a GPU is found
	try_run ( MIC_RUN_RESULT_VAR MIC_COMPILE_RESULT_VAR
	    ${CMAKE_BINARY_DIR}
	    ${CMAKE_MODULE_PATH}/Compiler/mic/testMICfound.cpp
	    CMAKE_FLAGS 
	    COMPILE_OUTPUT_VARIABLE MIC_COMPILE_OUTPUT_VAR
	    RUN_OUTPUT_VARIABLE MIC_RUN_OUTPUT_VAR )

	#message ( STATUS "MIC run output: ${MIC_COMPILE_OUTPUT_VAR} flag: ${MIC_COMPILE_RESULT_VAR}" )

	if    ( MIC_COMPILE_RESULT_VAR )

		message ( STATUS "MIC run output: ${MIC_RUN_OUTPUT_VAR} flag: ${MIC_RUN_RESULT_VAR}" )

		if     ( ${MIC_RUN_OUTPUT_VAR} MATCHES MIC )
			set ( USE_MIC TRUE )
		elseif ( ${MIC_RUN_OUTPUT_VAR} MATCHES HOST )
			set ( USE_MIC FALSE )
		else   ( )
			message ( FATAL_ERROR "Internal Error: MIC test examples return neighter MIC or HOST." )
		endif  ( )

	else  ( MIC_COMPILE_RESULT_VAR )
		# TODO: shouldn't be here an fatal error because testfile can't compile ?
		set ( USE_MIC FALSE )
	endif ( MIC_COMPILE_RESULT_VAR )

else  ( CMAKE_CXX_COMPILER_ID MATCHES Intel )
	set ( USE_MIC FALSE )
endif ( CMAKE_CXX_COMPILER_ID MATCHES Intel )