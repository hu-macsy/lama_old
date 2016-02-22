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

#include <scai/dmemo/Distribution.hpp>
#include <scai/dmemo/NoDistribution.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/test/TestMacros.hpp>

#include <scai/hmemo/HArray.hpp>
#include <scai/common/Math.hpp>

#include <boost/test/unit_test.hpp>

/** This routine is a rather general routine. It compares two arbitrary matrices
 *  (different distributions, different value types) whether the elements are
 *  close enough.
 *
 *  For efficiency on distributed matrices, it will do it row-wise.
 */

/**
 * @brief This method compares checks if two matrices have (nearly) the same values.
 *
 * @param[in] m1 is the first matrix
 * @param[in] m2 is the second matrix
 * @param[in] maximal absolute difference between two values
 * @param[in] tolerance in percentage that value might differ
 *
 * Both matrices must have the same dimensions, but they can have different types (double, float, ... )
 * or different storage formats ( DENSE, CSR, ELL, ... ) or different distributions.
 *
 * If the absolute difference between the values is greater than small than the values
 * must be in a tolerance that depends on the precision.
 */
static inline void testSameMatrix( const scai::lama::Matrix& m1, 
                                   const scai::lama::Matrix& m2, 
                                   scai::lama::Scalar small = scai::lama::Scalar( 0 ),
                                   scai::lama::Scalar tolerance = scai::lama::Scalar( 0.01 ) )
{
    using namespace scai;
    using namespace lama;

    const IndexType nRows = m1.getNumRows();
    const IndexType nCols = m1.getNumColumns();

    // require for same sizes, otherwise it is illegal to access elements

    BOOST_REQUIRE_EQUAL( nRows, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( nCols, m2.getNumColumns() );

    // create dense vectors for the rows with the same value type

    VectorCreateKeyType vectorType1( vectorformat::DENSE, m1.getValueType() );
    VectorCreateKeyType vectorType2( vectorformat::DENSE, m2.getValueType() );

    common::unique_ptr<Vector> ptrRow1( Vector::create( vectorType1 ) );
    common::unique_ptr<Vector> ptrRow2( Vector::create( vectorType2 ) );

    typedef double CompareType;  // complex type does not work

    CompareType tol = tolerance.getValue<CompareType>();

    // now compare all rows

    for ( IndexType i = 0; i < nRows; ++i )
    {
        // Note: rows will be broadcast in case of distributed matrices

        m1.getRow( *ptrRow1, i );
        m2.getRow( *ptrRow2, i );

        // compare the two vectors element-wise

        for ( IndexType j = 0; j < nCols; j++ )
        {
            // std::cout << i << ", " << j << ": " << ptrRow1->getValue( j ) << ", " << ptrRow2->getValue( j ) << std::endl;

            Scalar elem1 = ptrRow1->getValue( j );
            Scalar elem2 = ptrRow2->getValue( j );

            Scalar diff  = abs( elem1 - elem2 );

            if ( diff > small )
            {
                // for large numbers we just check for a tolerance

                BOOST_CHECK_CLOSE( elem1.getValue<CompareType>(), elem2.getValue<CompareType>(), tol );
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
void testSameMatrixStorage( const scai::lama::MatrixStorage<ValueType1>& m1, 
                            const scai::lama::MatrixStorage<ValueType2>& m2,
                            const scai::lama::Scalar small = 0.0,
                            const scai::lama::Scalar tolerance = 0.01 )
{
    const IndexType m = m1.getNumRows();
    const IndexType n = m1.getNumColumns();

    BOOST_REQUIRE_EQUAL( m, m2.getNumRows() );
    BOOST_REQUIRE_EQUAL( n, m2.getNumColumns() );

    scai::hmemo::HArray<ValueType1> row1(n);
    scai::hmemo::HArray<ValueType2> row2(n);

    typedef double CompareType;  // complex type does not work

    CompareType tolV   = tolerance.getValue<CompareType>();
    CompareType smallV = small.getValue<CompareType>();

    for ( IndexType i = 0; i < m; ++i )
    {
        m1.getRow( row1, i );
        m2.getRow( row2, i );

        // compare the two vectors element-wise

        scai::hmemo::ReadAccess<ValueType1> readRow1(row1);
        scai::hmemo::ReadAccess<ValueType2> readRow2(row2);

        for ( IndexType j = 0; j < n; j++ )
        {
            CompareType elem1 = static_cast<CompareType>( readRow1[j] );
            CompareType elem2 = static_cast<CompareType>( readRow2[j] );
       
            CompareType diff = scai::common::Math::abs( elem1 - elem2 );

            if ( diff >= smallV )
            {
                // if absolute difference is too big we check for close

                BOOST_CHECK_CLOSE( elem1, elem2, tolV );
            }
        }
    }
}
