###
 # @file scai/lamasolver/CMakeLists.txt
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
 # @brief Main specification file for CMake
 # @author Lauretta Schubert
 # @date 25.01.2016
 # @since 2.0.0
###

add_subdirectory( distributed )

scai_project ( 

	CLASSES

	HEADERS

	SOURCES

    	BiCGstabTest
        # BiCGTest  : disabled as long as MatrixTranspose for DenseMatrix is not supported
        CGTest
        # not yet: CGNETest
        CGNRTest
        # CGSTest: This test fails: testSolve, result not < 1E-4
        CommonLoggerTest
        CriterionTest
        FileLoggerTest
        GMRESTest
        InverseSolverTest
        IterationCountTest
        JacobiTest
        MINRESTest
        ResidualStagnationTest
        ResidualThresholdTest
        RichardsonTest
        SimpleAMGTest
        SolutionProxyTest
        SORTest
	    SpecializedJacobiTest
        TFQMRTest
        TimerTest
   )