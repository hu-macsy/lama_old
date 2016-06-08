/**
 * @file DIASparseMatrixTest.cpp
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
 * @brief Contains the implementation of the class DIASparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/matrix/DIASparseMatrix.hpp>
#include <scai/lama/test/SparseMatrixTest.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

extern bool base_test_case;
extern std::string testcase;

namespace scai
{
namespace lama
{
namespace DIASparseMatrixTest
{

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    DIASparseMatrix<ValueType> diaMatrix;
    SparseMatrixTest< DIASparseMatrix<ValueType> > diaSparseMatrixTest( diaMatrix );
    diaSparseMatrixTest.mMatrix.setContextPtr( loc );

    if ( base_test_case )
    {
        SPARSEMATRIX_COMMONTESTCASES( diaSparseMatrixTest );
    }
    else
    {
        diaSparseMatrixTest.runTests();
    }
}


template<typename ValueType>
void typeNameTest( )
{
    DIASparseMatrix<ValueType> diaMatrix;
    std::string s = diaMatrix.typeName();
    BOOST_CHECK( s.length() > 0 );
}

} /* end namespace DIASparseMatrixTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( DIASparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.DIASparseMatrixTest" )

LAMA_AUTO_TEST_CASE_CT( commonTestCases, DIASparseMatrixTest, scai::lama )
LAMA_AUTO_TEST_CASE_T( typeNameTest, DIASparseMatrixTest )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
