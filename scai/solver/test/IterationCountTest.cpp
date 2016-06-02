/**
 * @file IterationCountTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @brief Contains the implementation of the class IterationCountTest.
 * @author Alexander Büchel, Matthias Makulla
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/Jacobi.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

struct IterationCountTestConfig
{
    IterationCountTestConfig()
    {
        mCriterionDouble = new IterationCount( 5 );
        mCriterionFloat = new IterationCount( 10 );
    }

    ~IterationCountTestConfig()
    {
        delete mCriterionDouble;
        delete mCriterionFloat;
    }

    IterationCount* mCriterionDouble;
    IterationCount* mCriterionFloat;
};

BOOST_FIXTURE_TEST_SUITE( IterationCountTest, IterationCountTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.IterationCountTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getIterationExtrema(), 5 );
    BOOST_CHECK_EQUAL( mCriterionFloat->getIterationExtrema(), 10 );
    IterationCount* testcriterion = new IterationCount();
    BOOST_CHECK_EQUAL( testcriterion->getIterationExtrema(), 1 );
    IterationCount* testcriterion2 = new IterationCount( *testcriterion );
    BOOST_CHECK_EQUAL( testcriterion2->getIterationExtrema(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetIterationCountTest )
{
    BOOST_CHECK_EQUAL( 5, mCriterionDouble->getIterationExtrema() );
    mCriterionDouble->setIterationExtrema( 15 );
    BOOST_CHECK_EQUAL( 15, mCriterionDouble->getIterationExtrema() );
    BOOST_CHECK_EQUAL( 10, mCriterionFloat->getIterationExtrema() );
    mCriterionFloat->setIterationExtrema( 3 );
    BOOST_CHECK_EQUAL( 3, mCriterionFloat->getIterationExtrema() );
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
void testIsSatisfiedmethod( IndexType expectedIterationCount )
{
    Jacobi jacobi( "IterationCountTest Jacobi" );
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );
    DenseVector<ValueType> solution( 3, 1.0 );
    CriterionPtr criterion( new IterationCount( 10 ) );
    jacobi.initialize( coefficients );
    jacobi.setStoppingCriterion( criterion );
    jacobi.solve( solution, rhs );
    BOOST_CHECK_EQUAL( expectedIterationCount, jacobi.getIterationCount() );
}
;

BOOST_AUTO_TEST_CASE_TEMPLATE( testIsSatisfied, ValueType, scai_arithmetic_test_types )
{
    testIsSatisfiedmethod<ValueType>( 10 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    SCAI_COMMON_WRITEAT_TEST( *mCriterionDouble );
    SCAI_COMMON_WRITEAT_TEST( *mCriterionFloat );
    std::stringstream mStream;
    mStream << mCriterionDouble;
    std::string mString = mStream.str();
    BOOST_CHECK( mString.compare( "Maximal" ) );
    mStream.clear();
    mStream << mCriterionFloat;
    mString = mStream.str();
    BOOST_CHECK( mString.compare( "Minimal" ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
