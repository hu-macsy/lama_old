/**
 * @file NormTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Contains the implementation of the class NormTest
 * @author Alexander Büchel, Micha
 * @date 03.02.2012
 * @since 1.0.0
 */

#include <test/NormTest.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/expression/VectorExpressions.hpp>

LAMA_COMMON_TEST_CASE( NormTest, positiveHomogeneityTest )
	DenseVector<double> x( 4, 1.0 );
	Scalar s = 3.0;

	DenseVector<double> tmp( s* x );

	//Inequality test
	BOOST_CHECK_EQUAL( mNorm.apply ( tmp ), s* mNorm.apply( x ) );
LAMA_COMMON_TEST_CASE_END()

LAMA_COMMON_TEST_CASE( NormTest, triangleInequalityTest );
	DenseVector<double> x( 2, 2.0 );
	DenseVector<double> y( 2, 2.0 );
	DenseVector<double> z( x + y );

	BOOST_CHECK( mNorm.apply( z ) <= mNorm.apply( x ) + mNorm.apply( y ) );
LAMA_COMMON_TEST_CASE_END()

LAMA_COMMON_TEST_CASE( NormTest, ZeroVectorTest );
	DenseVector<double> x( 4, 0.0 );
	BOOST_CHECK_EQUAL( mNorm.apply( x ), 0.0 );
LAMA_COMMON_TEST_CASE_END()

LAMA_COMMON_TEST_CASE_RUNNER( NormTest )
{
    positiveHomogeneityTest();
    triangleInequalityTest();
    ZeroVectorTest();
}
