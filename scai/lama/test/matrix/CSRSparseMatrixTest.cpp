/**
 * @file test/matrix/CSRSparseMatrixTest.cpp
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
 * @brief Contains only test specific for the CSR Sparse matrix
 * @author Thomas Brandes
 * @date 24.03.2016
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/test/TestMacros.hpp>

#include <scai/common/TypeTraits.hpp>

using namespace scai;
using namespace lama;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CSRSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CSRSparseMatrixTest" );

/* ------------------------------------------------------------------------- */

typedef boost::mpl::list<ARITHMETIC_HOST> SCAI_ARITHMETIC_TYPES;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( defaultConstructorTest, ValueType, SCAI_ARITHMETIC_TYPES )
{
    CSRSparseMatrix<ValueType> matrix;

    // check zero sizes 

    BOOST_CHECK_EQUAL( 0, matrix.getNumRows() );
    BOOST_CHECK_EQUAL( 0, matrix.getNumColumns() );

    // check correct format / type

    BOOST_CHECK_EQUAL( common::TypeTraits<ValueType>::stype, matrix.getValueType() );
    BOOST_CHECK_EQUAL( Matrix::CSR, matrix.getFormat() );

    const CSRStorage<ValueType>& local = matrix.getLocalStorage();
    const CSRStorage<ValueType>& halo = matrix.getHaloStorage();

    BOOST_CHECK_EQUAL( local.getNumRows(), halo.getNumRows() );
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
