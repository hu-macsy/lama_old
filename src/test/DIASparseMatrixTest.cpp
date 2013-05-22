/**
 * @file DIASparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class DIASparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/matrix/DIASparseMatrix.hpp>

#include <test/SparseMatrixTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DIASparseMatrixTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.DIASparseMatrixTest" );

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, test_types ) {
    typedef T ValueType;

    DIASparseMatrix<ValueType> diaMatrix;
    SparseMatrixTest< DIASparseMatrix<ValueType> > diaSparseMatrixtest( diaMatrix );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in DIASparseMatrixTest." );
        SPARSEMATRIX_COMMONTESTCASES( diaSparseMatrixtest );
    }
    else
    {
        CONTEXTLOOP()
        {
            GETCONTEXT( context );
            diaSparseMatrixtest.mMatrix.setContext( context );
            LAMA_LOG_INFO( logger, "Using context = " << diaSparseMatrixtest.mMatrix.getContext().getType() );
            diaSparseMatrixtest.runTests();
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    DIASparseMatrix<double> diaMatrixd;
    std::string s = diaMatrixd.typeName();
    BOOST_CHECK_EQUAL( s, "DIASparseMatrix<double>" );

    DIASparseMatrix<float> diaMatrixf;
    s = diaMatrixf.typeName();
    BOOST_CHECK_EQUAL( s, "DIASparseMatrix<float>" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
