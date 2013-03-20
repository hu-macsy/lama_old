/**
 * @file ELLSparseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the implementation of the class ELLSparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 * $Id$
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/matrix/ELLSparseMatrix.hpp>

#include <test/SparseMatrixTest.hpp>
#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( ELLSparseMatrixTest )
;

LAMA_LOG_DEF_LOGGER( logger, "Test.ELLSparseMatrixTest" );

typedef boost::mpl::list<float,double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, T, test_types ) {
    typedef T ValueType;

    ELLSparseMatrix<ValueType> ellMatrix;
    SparseMatrixTest< ELLSparseMatrix<ValueType> > ellSparseMatrixtest( ellMatrix );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run method " << testcase << " in ELLSparseMatrixTest." );
        SPARSEMATRIX_COMMONTESTCASES( ellSparseMatrixtest );
    }
    else
    {
        CONTEXTLOOP()
        {
            GETCONTEXT( context );
            ellSparseMatrixtest.mMatrix.setContext( context );
            LAMA_LOG_INFO( logger, "Using context = " << ellSparseMatrixtest.mMatrix.getContext().getType() );
            ellSparseMatrixtest.runTests();
        }
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( typeNameTest )
{
    ELLSparseMatrix<double> ellMatrixd;
    std::string s = ellMatrixd.typeName();
    BOOST_CHECK_EQUAL( s, "ELLSparseMatrix<double>" );

    ELLSparseMatrix<float> ellMatrixf;
    s = ellMatrixf.typeName();
    BOOST_CHECK_EQUAL( s, "ELLSparseMatrix<float>" );
}
/* ------------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
