/**
 * @file ResidualStagnationTest.cpp
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
 * @brief Contains the implementation of the class ResidualStagnationTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/solver/criteria/ResidualStagnation.hpp>
#include <scai/solver/Jacobi.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/solver/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

struct ResidualStagnationTestConfig
{
    ResidualStagnationTestConfig()
    {
        NormPtr norm( new MaxNorm() );
        mCriterionDouble = new ResidualStagnation( norm, 6, Scalar( 10.0 ) );
        mCriterionFloat = new ResidualStagnation( norm, 4, Scalar( 10.0 ) );
    }

    ~ResidualStagnationTestConfig()
    {
        delete mCriterionDouble;
        delete mCriterionFloat;
    }

    ResidualStagnation* mCriterionDouble;
    ResidualStagnation* mCriterionFloat;
};

BOOST_FIXTURE_TEST_SUITE( ResidualStagnationTest, ResidualStagnationTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.ResidualStagnationTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getLookback(), 6 );
    BOOST_CHECK_EQUAL( mCriterionFloat->getLookback(), 4 );
    NormPtr norm( new MaxNorm() );
    ResidualStagnation* testcriterion = new ResidualStagnation( norm );
    BOOST_CHECK_EQUAL( testcriterion->getLookback(), 1 );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), 0.1 );
    ResidualStagnation* testcriterion2 = new ResidualStagnation( *testcriterion );
    BOOST_CHECK_EQUAL( testcriterion2->getLookback(), 1 );
    BOOST_CHECK_EQUAL( testcriterion2->getPrecision(), 0.1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    ResidualStagnation* testcriterion;
    testcriterion = ( ResidualStagnation* ) mCriterionFloat->copy();
    BOOST_CHECK_EQUAL( testcriterion->getLookback(), 4 );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), 10.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( NormPtrTest )
{
    NormPtr norm = mCriterionFloat->getNorm();
    //TODO: Check if this norm is really a L2Norm.
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetPrecisionTest )
{
    Scalar s = 6.0;
    NormPtr maxNorm( new MaxNorm() );
    ResidualStagnation* testcriterion = new ResidualStagnation( maxNorm );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), 0.1 );
    testcriterion->setPrecision( s );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), 6.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( GetAndSetLookBackTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getLookback(), 6 );
    BOOST_CHECK_EQUAL( mCriterionFloat->getLookback(), 4 );
    mCriterionDouble->setLookback( 4 );
    mCriterionFloat->setLookback( 2 );
    BOOST_CHECK_EQUAL( mCriterionDouble->getLookback(), 4 );
    BOOST_CHECK_EQUAL( mCriterionFloat->getLookback(), 2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( IsSatisfiedTest, ValueType, scai_arithmetic_test_types )
{
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get8x8SystemA<ValueType>();
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );
    DenseVector<ValueType> solution( 8, 0.0 );
    Jacobi solver( "StagnationTest solver" );
    NormPtr maxNorm( new MaxNorm() );
    CriterionPtr criterion( new ResidualStagnation( maxNorm, 2, Scalar( 1.1 ) ) );
    solver.setStoppingCriterion( criterion );
    solver.initialize( coefficients );
    solver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( 1, solver.getIterationCount() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    NormPtr maxNorm( new MaxNorm() );
    ResidualStagnation* testcriterion = new ResidualStagnation( maxNorm );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
