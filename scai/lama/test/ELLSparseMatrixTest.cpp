/**
 * @file ELLSparseMatrixTest.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class ELLSparseMatrixTest
 * @author Alexander Büchel
 * @date 02.03.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/lama/matrix/ELLSparseMatrix.hpp>
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
namespace ELLSparseMatrixTest
{

template<typename ValueType>
void commonTestCases( ContextPtr loc )
{
    ELLSparseMatrix<ValueType> ellMatrix;
    SparseMatrixTest< ELLSparseMatrix<ValueType> > ellSparseMatrixTest( ellMatrix );
    ellSparseMatrixTest.mMatrix.setContextPtr( loc );

    if ( base_test_case )
    {
        SPARSEMATRIX_COMMONTESTCASES( ellSparseMatrixTest );
    }
    else
    {
        ellSparseMatrixTest.runTests();
    }
}

template<typename ValueType>
void typeNameTest( )
{
    ELLSparseMatrix<ValueType> ellMatrix;
    std::string s = ellMatrix.typeName();
    BOOST_CHECK( s.length() > 0 );
}

} /* end namespace ELLSparseMatrixTest */

} /* end namespace lama */

} /* end namespace scai */

/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE( ELLSparseMatrixTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.ELLSparseMatrixTest" )

LAMA_AUTO_TEST_CASE_CT( commonTestCases, ELLSparseMatrixTest, scai::lama )
LAMA_AUTO_TEST_CASE_T( typeNameTest, ELLSparseMatrixTest )
/* ------------------------------------------------------------------------------------------------------------------ */

BOOST_AUTO_TEST_SUITE_END()
