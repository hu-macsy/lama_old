/**
 * @file P_COOSparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class P_COOSparseMatrixTest
 * @author Alexander Büchel
 * @date 10.05.2012
 * @since 1.0.0
 */

#include <test/distributed/P_SparseMatrixTest.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/COOSparseMatrix.hpp>
#include <test/TestMacros.hpp>

using namespace lama;
using namespace memory;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( P_COOSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest.P_COOSparseMatrixTest" );

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, ValueType, test_types )
{
    P_SparseMatrixTest<COOSparseMatrix<ValueType> > p_cooSparseMatrixtest;

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run test method " << testcase << " in P_COOSparseMatrixTest." );
        PSPARSEMATRIXTEST_COMMONTESTCASES( p_cooSparseMatrixtest );
    }
    else
    {
        p_cooSparseMatrixtest.runTests();
    }
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
