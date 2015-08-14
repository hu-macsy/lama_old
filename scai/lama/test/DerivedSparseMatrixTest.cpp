/**
 * @file CSRSparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class CSRSparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/matrix/ELLSparseMatrix.hpp>
#include <scai/lama/matrix/JDSSparseMatrix.hpp>

#include <test/SparseMatrixTest.hpp>

using namespace lama;
using namespace memory;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DerivedSparseMatrixTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest.DerivedSparseMatrixTest" )

typedef boost::mpl::list<CSRSparseMatrix<float>, CSRSparseMatrix<double>, COOSparseMatrix<float>, COOSparseMatrix<double>,
        DIASparseMatrix<float>, DIASparseMatrix<double>, ELLSparseMatrix<float>, ELLSparseMatrix<double>,
        JDSSparseMatrix<float>, JDSSparseMatrix<double> > test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( clear, MatrixType, test_types )
{
    MatrixType matrix;
    SparseMatrixTest<MatrixType> sparseMatrixtest( matrix );
    sparseMatrixtest.clearTest();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( conversions, MatrixType, test_types )
{
    MatrixType matrix;
    SparseMatrixTest<MatrixType> sparseMatrixtest( matrix );
    sparseMatrixtest.testConversions();
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
