/**
 * @file SameMatrixHelper.cpp
 *
 * @license
 * Copyright (c) 2012
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
 * @brief Implementaton of methods to verify that two matrices have same values.
 * @author Thomas Brandes
 * @date 17.06.2012
 * $Id$
 */

#include <test/SameMatrixHelper.hpp>
#include <test/TestMacros.hpp>

using namespace lama;

/************************************************************************
 *  Planned version (avoids use of templates)                            *
 ************************************************************************/

void verifySameMatrixAll( const Matrix& m1, const Matrix& m2 )
{
    const IndexType n = m1.getNumRows();

    BOOST_CHECK_EQUAL( n, m1.getNumColumns() );
    BOOST_CHECK_EQUAL( n, m2.getNumRows() );
    BOOST_CHECK_EQUAL( n, m2.getNumColumns() );

    VectorPtr ptrX1 = m1.createDenseVector( m1.getColDistributionPtr(), 1.0 );
    VectorPtr ptrX2 = m2.createDenseVector( m2.getColDistributionPtr(), 1.0 );

    VectorPtr ptrY1 = m1.createDenseVector( m1.getDistributionPtr(), 0.0 );
    VectorPtr ptrY2 = m2.createDenseVector( m2.getDistributionPtr(), 0.0 );

//    Vector& x1 = *ptrX1;
//    Vector& x2 = *ptrX2;
    Vector& y1 = *ptrY1;
    Vector& y2 = *ptrY2;

    // This does not work
    // y1 = m1 * x1;
    // y2 = m2 * x2;

    for( IndexType i = 0; i < n; i++ )
    {
        Scalar s1 = y1.getValue( i );
        Scalar s2 = y2.getValue( i );
        BOOST_CHECK_CLOSE( s1.getValue<double>(), s2.getValue<double>(), 1.0e-5 );
    }
}

/************************************************************************
 *  Planned version  (works on columns or rows)                          *
 ************************************************************************/

void assertSameMatrix( const Matrix& m1, const Matrix& m2 )
{
    std::cout << "compare m1 = " << m1 << " with m2 = " << m2 << std::endl;

    const IndexType m = m1.getNumRows();
    const IndexType n = m1.getNumColumns();

    // require for same sizes, otherwise it is illegal to access elements

    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    DistributionPtr replicated( DistributionPtr( new NoDistribution( m ) ) );

    // create replicated vectors for the rows with the same value type

    VectorPtr ptrRow1 = m1.createDenseVector( replicated, 0.0 );
    VectorPtr ptrRow2 = m2.createDenseVector( replicated, 0.0 );

    // now compare all rows

    for( IndexType i = 0; i < m; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices

        std::cout << "compare row " << i << std::endl;

        m1.getRow( *ptrRow1, i );
        m2.getRow( *ptrRow2, i );

        std::cout << "compare row element-wise: " << i << std::endl;

        // compare the two vectors element-wise

        for( IndexType j = 0; j < n; j++ )
        {
            Scalar s1 = ptrRow1->getValue( j );
            Scalar s2 = ptrRow2->getValue( j );
            LAMA_CHECK_SCALAR_CLOSE( s1, s2, double, 1 );
        }
    }
}
