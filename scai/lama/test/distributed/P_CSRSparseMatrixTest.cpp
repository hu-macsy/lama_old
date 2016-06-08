/**
 * @file P_CSRSparseMatrixTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Affero General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Contains the implementation of the class P_CSRSparseMatrixTest
 * @author Alexander Büchel
 * @date 10.05.2012
 */

#include <scai/lama/test/distributed/P_SparseMatrixTest.hpp>

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( P_CSRSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SparseMatrixTest.P_CSRSparseMatrixTest" );

typedef boost::mpl::list<float, double> test_types;

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( commonTestCases, ValueType, test_types )
{
    P_SparseMatrixTest< CSRSparseMatrix<ValueType> > p_csrSparseMatrixtest;

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run test method " << testcase << " in P_CSRSparseMatrixTest." );
        PSPARSEMATRIXTEST_COMMONTESTCASES( p_csrSparseMatrixtest );
    }
    else
    {
        p_csrSparseMatrixtest.runTests();
    }
}

/* ------------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
