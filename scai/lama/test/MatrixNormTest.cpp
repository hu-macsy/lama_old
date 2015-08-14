/**
 * @file MatrixNormTest.cpp
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
 * @author Eric Schricker
 * @date 29.04.2015
 * @since 1.1.0
 */

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/Scalar.hpp>
#include <scai/lama/matrix/Matrix.hpp>
#include <scai/lama/matrix/all.hpp>
#include <scai/lama/norm/Norm.hpp>

#include <test/TestMacros.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

typedef boost::mpl::list<float, double> test_types;

using namespace lama;
using namespace memory;

BOOST_AUTO_TEST_SUITE( MatrixNormTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MatrixNormTest" )

template<typename MatrixType>
void l1NormTestMethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

	ValueType tmp = 16;

	Scalar result(tmp);

	MatrixType matrix;
	matrix.setIdentity(8);

	matrix.setContext( context );

	matrix *= 2.0;

	Scalar l1Norm = matrix.l1Norm();

	BOOST_CHECK_EQUAL( l1Norm, result );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( l1NormTest, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        l1NormTestMethod< CSRSparseMatrix<ValueType> >( context );
        l1NormTestMethod< ELLSparseMatrix<ValueType> >( context );
        l1NormTestMethod< COOSparseMatrix<ValueType> >( context );
        l1NormTestMethod< JDSSparseMatrix<ValueType> >( context );
        l1NormTestMethod< DIASparseMatrix<ValueType> >( context );
        l1NormTestMethod< DenseMatrix<ValueType> >( context );
        l1NormTestMethod< DIASparseMatrix<ValueType> >( context );
    }
}

template<typename MatrixType>
void l2NormTestMethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

	ValueType tmp = 32.0;

	Scalar result = sqrt(tmp);

	MatrixType matrix;
	matrix.setIdentity(8);

	matrix.setContext( context );

	matrix *= 2.0;

	Scalar l2Norm = matrix.l2Norm();

	BOOST_CHECK_EQUAL( l2Norm, result );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( l2NormTest, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        l2NormTestMethod< CSRSparseMatrix<ValueType> >( context );
        l2NormTestMethod< ELLSparseMatrix<ValueType> >( context );
        l2NormTestMethod< COOSparseMatrix<ValueType> >( context );
        l2NormTestMethod< JDSSparseMatrix<ValueType> >( context );
        l2NormTestMethod< DIASparseMatrix<ValueType> >( context );
        l2NormTestMethod< DenseMatrix<ValueType> >( context );
        l2NormTestMethod< DIASparseMatrix<ValueType> >( context );
    }
}

template<typename MatrixType>
void maxNormTestMethod( ContextPtr context )
{
    typedef typename MatrixType::MatrixValueType ValueType;

	ValueType tmp = 2;

	Scalar result(tmp);

	MatrixType matrix;
	matrix.setIdentity(8);

	matrix.setContext( context );

	matrix *= 2.0;

	Scalar maxNorm = matrix.maxNorm();

	BOOST_CHECK_EQUAL( maxNorm, result );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( maxNormTest, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );
        maxNormTestMethod< CSRSparseMatrix<ValueType> >( context );
        maxNormTestMethod< ELLSparseMatrix<ValueType> >( context );
        maxNormTestMethod< COOSparseMatrix<ValueType> >( context );
        maxNormTestMethod< JDSSparseMatrix<ValueType> >( context );
        maxNormTestMethod< DIASparseMatrix<ValueType> >( context );
        maxNormTestMethod< DenseMatrix<ValueType> >( context );
        maxNormTestMethod< DIASparseMatrix<ValueType> >( context );
    }
}

BOOST_AUTO_TEST_SUITE_END();
