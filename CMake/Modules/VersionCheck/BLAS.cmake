###
 # @file BLAS.cmake
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
 # @brief Check BLAS Version
 # @author Jan Ecker
 # @date 25.04.2013
 # @since 1.0.0
###

### MKL

if    ( SCAI_BLAS_NAME MATCHES MKL )
	try_run ( MKL_RUN_RESULT_VAR MKL_COMPILE_RESULT_VAR
		    ${CMAKE_BINARY_DIR}
	    	${CMAKE_MODULE_PATH}/VersionCheck/mkl.cpp
	    	CMAKE_FLAGS 
	    	-DINCLUDE_DIRECTORIES:STRING=${SCAI_SCAI_BLAS_INCLUDE_DIR}
	    	COMPILE_OUTPUT_VARIABLE MKL_COMPILE_OUTPUT_VAR
	    	RUN_OUTPUT_VARIABLE MKL_RUN_OUTPUT_VAR )

    set ( BLAS_VERSION ${MKL_RUN_OUTPUT_VAR} )
endif ( SCAI_BLAS_NAME MATCHES MKL )

### BLAS

if    ( SCAI_BLAS_NAME MATCHES BLAS )
	get_filename_component ( _LIB_PATH ${BLAS_blas_LIBRARY} PATH )
	file ( GLOB _LIBRARIES ${_LIB_PATH}/libblas.so.* )

	## look for lib with longest name (version number)
	set ( LEN 0 )
	foreach    ( _LIB ${_LIBRARIES} )
		string ( LENGTH ${_LIB} TMP_LEN )
		if    ( ${TMP_LEN} GREATER ${LEN} )
			set ( LEN ${TMP_LEN} )
			set ( VAL ${_LIB} )
		endif ( ${TMP_LEN} GREATER ${LEN} )
	endforeach ( _LIB ${_LIBRARIES} )

	if    ( DEFINED VAL )
    	string ( REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+)" BLAS_VERSION ${VAL} )
    else  ( DEFINED VAL )
    	set ( BLAS_VERSION "unknown BLAS version")
    endif ( DEFINED VAL )
endif ( SCAI_BLAS_NAME MATCHES BLAS )

### INTERNAL BLAS

if    ( SCAI_BLAS_NAME MATCHES INTERNALBLAS )
    set ( BLAS_VERSION ${SCAI_BLASKERNEL_VERSION} )
endif ( SCAI_BLAS_NAME MATCHES INTERNALBLAS )

