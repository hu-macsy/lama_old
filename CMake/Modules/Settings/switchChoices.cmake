###
 # @file switchChoices.cmake
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
 # @brief Sets all lists of choices that can be switch in the LAMA configutation and defines their default behaviour.
 # @author Lauretta Schubert
 # @date 01.04.2016
 # @since 2.0.0
###

set ( TRUE_FALSE_CHOICES ON OFF )

set ( BUILD_DOC_DEFAULT ON )
set ( BUILD_EXAMPLES_DEFAULT ON )
set ( BUILD_TEST_DEFAULT ON )
set ( USE_CODE_COVERAGE_DEFAULT OFF )

set ( CMAKE_BUILD_TYPE_CHOICES "None" "Debug" "Release" "RelWithDebInfo" "MinSizeRel" ) 
set ( CMAKE_BUILD_TYPE_DEFAULT "Debug" )

set ( SCAI_ASSERT_CHOICES DEBUG ERROR OFF )
set ( SCAI_ASSERT_DEFAULT DEBUG )

set ( SCAI_BLAS_LIBRARY_CHOICES auto MKL BLAS INTERNALBLAS )
set ( SCAI_BLAS_LIBRARY_DEFAULT auto )

set ( SCAI_DOC_TYPE_CHOICES html json xml )
set ( SCAI_DOC_TYPE_DEFAULT html )

set ( SCAI_LIBRARY_TYPE_CHOICES STATIC SHARED )
set ( SCAI_LIBRARY_TYPE_DEFAULT SHARED )

set ( SCAI_LOGGING_CHOICES TRACE DEBUG INFO WARN ERROR OFF )
# no default, decision depending on choosen CMAKE_BUILD_TYPE --> see Settings/logging
