/**
 * @file JDSSparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class JDSSparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/matrix/JDSSparseMatrix.hpp>

#include <test/SparseMatrixTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( JDSSparseMatrixTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest.JDSSparseMatrixTest" )

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, test_types )
{
    typedef T ValueType;

    JDSSparseMatrix<ValueType> jdsMatrix;
    SparseMatrixTest< JDSSparseMatrix<ValueType> > jdsSparseMatrixtest( jdsMatrix );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in JDSSparseMatrixTest." );
        SPARSEMATRIX_COMMONTESTCASES( jdsSparseMatrixtest );
    }
    else
    {
        CONTEXTLOOP()
        {
            GETCONTEXT( context );
            jdsSparseMatrixtest.mMatrix.setContext( context );
            LAMA_LOG_INFO( logger, "Using context = " << jdsSparseMatrixtest.mMatrix.getContext().getType() );
            jdsSparseMatrixtest.runTests();
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    JDSSparseMatrix<double> jdsMatrixd;
    std::string s = jdsMatrixd.typeName();
    BOOST_CHECK_EQUAL( s, "JDSSparseMatrix<double>" );

    JDSSparseMatrix<float> jdsMatrixf;
    s = jdsMatrixf.typeName();
    BOOST_CHECK_EQUAL( s, "JDSSparseMatrix<float>" );
}
/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
