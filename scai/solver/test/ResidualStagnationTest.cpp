/**
 * @file ResidualStagnationTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Test routines for the class ResidualStagnation.
 * @author Thomas Brandes, Malte Foerster
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

using namespace scai;
using namespace solver;
using namespace lama;
using namespace hmemo;

/* --------------------------------------------------------------------- */

typedef SCAI_TEST_TYPE ValueType;

struct ResidualStagnationTestConfig
{
    ResidualStagnationTestConfig()
    {
        NormPtr<ValueType> norm( new MaxNorm<ValueType>() );
        mCriterion1 = new ResidualStagnation<ValueType>( norm, 6, ValueType( 10 ) );
        mCriterion2 = new ResidualStagnation<ValueType>( norm, 4, ValueType( 10 ) );
    }

    ~ResidualStagnationTestConfig()
    {
        delete mCriterion1;
        delete mCriterion2;
    }

    ResidualStagnation<ValueType>* mCriterion1;
    ResidualStagnation<ValueType>* mCriterion2;
};

BOOST_FIXTURE_TEST_SUITE( ResidualStagnationTest, ResidualStagnationTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.ResidualStagnationTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mCriterion1->getLookback(), IndexType( 6 ) );
    BOOST_CHECK_EQUAL( mCriterion2->getLookback(), IndexType( 4 ) );
    NormPtr<ValueType> norm( new MaxNorm<ValueType>() );
    ResidualStagnation<ValueType>* testcriterion = new ResidualStagnation<ValueType>( norm );
    BOOST_CHECK_EQUAL( testcriterion->getLookback(), IndexType( 1 ) );
    ValueType epsDefault = 0.1;
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), epsDefault );
    ResidualStagnation<ValueType>* testcriterion2 = new ResidualStagnation<ValueType>( *testcriterion );
    BOOST_CHECK_EQUAL( testcriterion2->getLookback(), IndexType( 1 ) );
    BOOST_CHECK_EQUAL( testcriterion2->getPrecision(), epsDefault );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    std::unique_ptr<ResidualStagnation<ValueType> > testcriterion( mCriterion2->copy() );
    BOOST_CHECK_EQUAL( testcriterion->getLookback(), IndexType( 4 ) );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), 10.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( NormPtrTest )
{
    NormPtr<ValueType> norm = mCriterion2->getNorm();
    //TODO: Check if this norm is really a L2Norm.
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetPrecisionTest )
{
    ValueType s = 6;
    NormPtr<ValueType> maxNorm( new MaxNorm<ValueType>() );
    ResidualStagnation<ValueType>* testcriterion = new ResidualStagnation<ValueType>( maxNorm );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), ValueType( 0.1 ) );
    testcriterion->setPrecision( s );
    BOOST_CHECK_EQUAL( testcriterion->getPrecision(), s );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( GetAndSetLookBackTest )
{
    BOOST_CHECK_EQUAL( mCriterion1->getLookback(), IndexType( 6 ) );
    BOOST_CHECK_EQUAL( mCriterion2->getLookback(), IndexType( 4 ) );
    mCriterion1->setLookback( 4 );
    mCriterion2->setLookback( 2 );
    BOOST_CHECK_EQUAL( mCriterion1->getLookback(), IndexType( 4 ) );
    BOOST_CHECK_EQUAL( mCriterion2->getLookback(), IndexType( 2 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( IsSatisfiedTest, ValueType, scai_numeric_test_types )
{
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get8x8SystemA<ValueType>();
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );
    auto solution = denseVectorFill<ValueType>( 8, 0 );
    Jacobi<ValueType> solver( "StagnationTest solver" );
    NormPtr<ValueType> maxNorm( new MaxNorm<ValueType>() );
    CriterionPtr<ValueType> criterion( new ResidualStagnation<ValueType>( maxNorm, 2, ValueType( 1.1 ) ) );
    solver.setStoppingCriterion( criterion );
    solver.initialize( coefficients );
    solver.solve( solution, rhs );
    BOOST_CHECK_EQUAL( IndexType( 1 ), solver.getIterationCount() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    NormPtr<ValueType> maxNorm( new MaxNorm<ValueType>() );
    ResidualStagnation<ValueType>* testcriterion = new ResidualStagnation<ValueType>( maxNorm );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
