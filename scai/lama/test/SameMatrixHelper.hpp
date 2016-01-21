/**
 * @file SameMatrixHelper.hpp
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
 * @brief This file contains routines to check the equality of distributed matrices
 * @author Thomas Brandes
 * @date 18.06.2012
 * @since 1.0.0
 */

#pragma once

#include <scai/lama/matrix/Matrix.hpp>

#include <scai/lama/distribution/Distribution.hpp>
#include <scai/lama/distribution/NoDistribution.hpp>
#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <scai/hmemo/HArray.hpp>

#include <boost/test/unit_test.hpp>

/** This routine is a rather general routine. It compares two arbitrary matrices
 *  (different distributions, different value types) whether the elements are
 *  close enough.
 *
 *  For efficiency on distributed matrices, it will do it row-wise.
 */

using namespace scai::lama;
using namespace scai::hmemo;


/**
 * @brief testSameMatrix() checks whether two matrices have the same dimensions and values. Equality of values is checked
 * by SCAI_CHECK_CLOSE, zero values have to be zero in both matrices (check testSameMatrixClose as alternative)
 */
template<typename MatrixType1, typename MatrixType2>
void testSameMatrix( const MatrixType1& m1, const MatrixType2& m2 )
{
    typedef typename MatrixType1::MatrixValueType ValueType1;
    typedef typename MatrixType2::MatrixValueType ValueType2;

    const IndexType m = m1.getNumRows();
    const IndexType n = m1.getNumColumns();
    // require for same sizes, otherwise it is illegal to access elements
    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    // create replicated vectors for the rows with the same value type
    VectorPtr ptrRow1(DenseVector<ValueType1>::createVector(m1.getValueType(), DistributionPtr( new NoDistribution( m ))));
    VectorPtr ptrRow2(DenseVector<ValueType2>::createVector(m2.getValueType(), DistributionPtr( new NoDistribution( m ))));

    // now compare all rows
    for ( IndexType i = 0; i < m; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices
        m1.getRow( *ptrRow1, i );
        m2.getRow( *ptrRow2, i );
        // compare the two vectors element-wise

        //#pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
        for ( IndexType j = 0; j < n; j++ )
        {
            SCAI_CHECK_CLOSE( ptrRow1->getValue( j ), ptrRow2->getValue( j ), eps<ValueType1>() )
        }
    }
}

/**
 * @brief testSameMatrixClose() checks whether two matrices have the same dimensions and values. Equality of values is
 * checked by SCAI_CHECK_CLOSE, if a value in the expected matrix is small (or zero) SCAI_CHECK_SMALL is used instead
 */
template<typename MatrixType1, typename MatrixType2>
void testSameMatrixClose( const MatrixType1& m1, const MatrixType2& m2 )
{
    typedef typename MatrixType1::MatrixValueType ValueType1;
    typedef typename MatrixType2::MatrixValueType ValueType2;

    const IndexType m = m1.getNumRows();
    const IndexType n = m1.getNumColumns();
    // require for same sizes, otherwise it is illegal to access elements
    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    // create replicated vectors for the rows with the same value type
    VectorPtr ptrRow1(DenseVector<ValueType1>::createVector(m1.getValueType(), DistributionPtr( new NoDistribution( m ))));
    VectorPtr ptrRow2(DenseVector<ValueType2>::createVector(m2.getValueType(), DistributionPtr( new NoDistribution( m ))));

    // now compare all rows
    for ( IndexType i = 0; i < m; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices
        m1.getRow( *ptrRow1, i );
        m2.getRow( *ptrRow2, i );
        // compare the two vectors element-wise

        //#pragma omp parallel for schedule(SCAI_OMP_SCHEDULE)
        for ( IndexType j = 0; j < n; j++ )
        {
            if ( abs( ptrRow1->getValue( j ) ) < small<ValueType1>() )
            {
                SCAI_CHECK_SCALAR_SMALL( ptrRow2->getValue( j ), ValueType2, small<ValueType2>() )
            }
            else
            {
                SCAI_CHECK_CLOSE( ptrRow1->getValue( j ), ptrRow2->getValue( j ), eps<ValueType1>() )
            }
        }
    }
}

/**
 * @brief testSameMatrixStorage() checks whether two matrices have the same dimensions and values. Equality of values is
 * checked by SCAI_CHECK_CLOSE, zero values have to be zero in both matrices (check testSameMatrixStorageClose as
 * alternative)
 */
template<typename ValueType1, typename ValueType2>
void testSameMatrixStorage( const MatrixStorage<ValueType1>& m1, const MatrixStorage<ValueType2>& m2 )
{
    const IndexType m = m1.getNumRows();
    const IndexType n = m1.getNumColumns();

    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    HArray<ValueType1> row1(n);
    HArray<ValueType2> row2(n);

    for ( IndexType i = 0; i < m; ++i )
    {
        m1.getRow( row1, i );
        m2.getRow( row2, i );
        // compare the two vectors element-wise

        ReadAccess<ValueType1> readRow1(row1);
        ReadAccess<ValueType2> readRow2(row2);
        for ( IndexType j = 0; j < n; j++ )
        {
            SCAI_CHECK_CLOSE( (readRow1.get())[j], (readRow2.get())[j], eps<ValueType1>() )
        }
    }
}

/**
 * @brief testSameMatrixStorageClose() checks whether two matrices have the same dimensions and values. Equality of
 * values is checked by SCAI_CHECK_CLOSE, if a value in the expected matrix is small (or zero) SCAI_CHECK_SMALL is used
 * instead.
 */
template<typename ValueType1, typename ValueType2>
void testSameMatrixStorageClose( const MatrixStorage<ValueType1>& m1, const MatrixStorage<ValueType2>& m2 )
{
    const IndexType m = m1.getNumRows();
    const IndexType n = m1.getNumColumns();

    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    HArray<ValueType1> row1(n);
    HArray<ValueType2> row2(n);

    for ( IndexType i = 0; i < m; ++i )
    {
        m1.getRow( row1, i );
        m2.getRow( row2, i );
        // compare the two vectors element-wise

        ReadAccess<ValueType1> readRow1(row1);
        ReadAccess<ValueType2> readRow2(row2);
        for ( IndexType j = 0; j < n; j++ )
        {
            if ( abs( (readRow1.get())[j] ) < small<ValueType1>() )
            {
                SCAI_CHECK_SCALAR_SMALL( Scalar((readRow2.get())[j]), ValueType2, small<ValueType2>() )
            }
            else
            {
                SCAI_CHECK_CLOSE( (readRow1.get())[j], (readRow2.get())[j], eps<ValueType1>() )
            }
        }
    }
}

