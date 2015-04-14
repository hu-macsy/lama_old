/**
 * @file SameMatrixHelper.hpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief This file contains routines to check the equality of distributed matrices
 * @author Thomas Brandes
 * @date 18.06.2012
 * @since 1.0.0
 */
#ifndef LAMA_SAME_MATRIX_HELPER_HPP_
#define LAMA_SAME_MATRIX_HELPER_HPP_

#include <lama/matrix/Matrix.hpp>

#include <lama/distribution/Distribution.hpp>
#include <lama/distribution/NoDistribution.hpp>
#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>

#include <lama/DenseVector.hpp>

#include <boost/test/unit_test.hpp>

/** This routine is a rather general routine. It compares two arbitrary matrices
 *  (different distributions, different value types) whether the elements are
 *  close enough.
 *
 *  The comparison is done by multiplying each matrix with a unity vector.
 */

template<typename MatrixType1,typename MatrixType2>
void testSameMatrix( const MatrixType1& m1, const MatrixType2& m2 )
{
    typedef typename MatrixType1::MatrixValueType ValueType1;
    typedef typename MatrixType2::MatrixValueType ValueType2;

    // Both matrices must be matrices of the same size

    const lama::IndexType m = m1.getNumRows();
    const lama::IndexType n = m1.getNumColumns();

    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    lama::DenseVector<ValueType1> x1( m1.getColDistributionPtr(), 1.0 );
    lama::DenseVector<ValueType2> x2( m2.getColDistributionPtr(), 1.0 );

    lama::DenseVector<ValueType1> y1( m1.getDistributionPtr(), 0.0 );
    lama::DenseVector<ValueType2> y2( m2.getDistributionPtr(), 0.0 );

    y1 = m1 * x1; // YA = mA * XA;
    y2 = m2 * x2; // YB = mB * XB;

    // replicate the vectors to avoid communication overhead for single elements

    lama::DistributionPtr replicated( new lama::NoDistribution( m ) );

    y1.redistribute( replicated );
    y2.redistribute( replicated );

    for ( lama::IndexType i = 0; i < m; i++ )
    {
        lama::Scalar s1 = y1.getValue( i );
        lama::Scalar s2 = y2.getValue( i );
        BOOST_CHECK_CLOSE( s1.getValue<double>(), s2.getValue<double>(), 0.01 );
    }
}

/** This test compares two general matrices by comparing all matrix elements.
 *
 *  For efficiency on distributed matrices, it will do it row-wise.
 */

void assertSameMatrix( const lama::Matrix& m1, const lama::Matrix& m2 );

#endif // LAMA_SAME_MATRIX_HELPER_HPP_
