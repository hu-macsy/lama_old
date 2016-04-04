###
 # @file scai/tracing/test/Makefile
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
 # @brief Makefile for building the executables required by the test script
 # @author Jan Ecker
 # @date 18.02.2016
 # @since 2.0.0
###

include make.inc

INCLUDE_FLAGS = ${SCAI_INCLUDES} -I${SCAI_BOOST_INCLUDE_DIR}

CXX_FLAGS = ${SCAI_CXX_FLAGS}

LD_FLAGS = ${SCAI_LD_FLAGS}

DEFINES =


simple: simpleTracing.exe

.SUFFIXES: .cpp .o
.SUFFIXES: .o .exe

.o.exe:
	${CXX} -o $*.exe $*.o ${LD_FLAGS}

.cpp.o:
	${CXX} ${CXX_FLAGS} ${DEFINES} ${INCLUDE_FLAGS} -c $*.cpp

clean:
	rm -rf *.exe *.o *.ct.* *.time
