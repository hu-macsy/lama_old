/**
 * @file ResidualStagnationTest.cpp
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
 * @brief Contains the implementation of the class ResidualStagnationTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 14.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/solver/criteria/ResidualStagnation.hpp>
#include <scai/lama/solver/DefaultJacobi.hpp>
#include <scai/lama/norm/MaxNorm.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace lama;
using namespace memory;

typedef boost::mpl::list<double, float> test_types;

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

BOOST_AUTO_TEST_CASE_TEMPLATE( IsSatisfiedTest, ValueType, test_types )
{
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get8x8SystemA<ValueType>();
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );
    DenseVector<ValueType> solution( 8, 0.0 );
    DefaultJacobi solver( "StagnationTest solver" );
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
    LAMA_WRITEAT_PTR_TEST( testcriterion );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
