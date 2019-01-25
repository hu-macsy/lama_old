/**
 * @file ResidualThresholdTest.cpp
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
 * @brief Test routines for the class ResidualThreshold.
 * @author Malte Foerster
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>

#include <scai/lama/matrix/CSRSparseMatrix.hpp>

#include <scai/lama/norm/MaxNorm.hpp>

#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/Jacobi.hpp>
#include <scai/lama/norm/L2Norm.hpp>

#include <scai/lama/expression/MatrixVectorExpressions.hpp>
#include <scai/lama/expression/VectorExpressions.hpp>

#include <scai/solver/test/EquationHelper.hpp>
#include <scai/solver/test/TestMacros.hpp>

#include <memory>

using namespace scai;
using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

typedef SCAI_TEST_TYPE ValueType;

/* --------------------------------------------------------------------- */

struct ResidualThresholdTestConfig
{
    ResidualThresholdTestConfig()
    {
        NormPtr<ValueType> mNorm( new L2Norm<ValueType>() );
        mCriterionDouble.reset( new ResidualThreshold<ValueType>( mNorm, ValueType( 1e-6 ), ResidualCheck::Absolute ) );
        mCriterionFloat.reset( new ResidualThreshold<ValueType>( mNorm, ValueType( 1e-4 ), ResidualCheck::Relative ) );
    }

    ~ResidualThresholdTestConfig()
    {
    }

    std::unique_ptr<ResidualThreshold<ValueType> > mCriterionDouble;
    std::unique_ptr<ResidualThreshold<ValueType> > mCriterionFloat;
};

BOOST_FIXTURE_TEST_SUITE( ResidualThresholdTest, ResidualThresholdTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.ResidualThresholdTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getPrecision(), ValueType( 1e-6 ) );
    BOOST_CHECK_EQUAL( mCriterionDouble->getCheckMode(), ResidualCheck::Absolute );
    BOOST_CHECK_EQUAL( mCriterionFloat->getPrecision(), ValueType( 1e-4 ) );
    BOOST_CHECK_EQUAL( mCriterionFloat->getCheckMode(), ResidualCheck::Relative );

    // testing for the default settings

    ValueType defaultEps = 1e-5;
    const ResidualThreshold<ValueType>* testcriterion1 = new ResidualThreshold<ValueType>();
    BOOST_CHECK_EQUAL( testcriterion1->getPrecision(), defaultEps );
    BOOST_CHECK_EQUAL( testcriterion1->getCheckMode(), ResidualCheck::Relative );
    BOOST_CHECK_EQUAL( testcriterion1->getFirstNormResult(), ValueType( -1 ) );

    NormPtr<ValueType> testnorm( new MaxNorm<ValueType>() );
    std::unique_ptr<ResidualThreshold<ValueType> > testcriterion2( new ResidualThreshold<ValueType>( testnorm ) );
    BOOST_CHECK_EQUAL( testcriterion2->getFirstNormResult(), -1.0 );
    BOOST_CHECK_EQUAL( testcriterion2->getPrecision(), ValueType( 1e-5 ) );
    BOOST_CHECK_EQUAL( testcriterion2->getCheckMode(), ResidualCheck::Relative );

    // verify copy constructor

    std::unique_ptr<ResidualThreshold<ValueType> > testcriterion3( new ResidualThreshold<ValueType>( *testcriterion1 ) );
    BOOST_CHECK_EQUAL( testcriterion3->getFirstNormResult(), ValueType( -1 ) );
    BOOST_CHECK_EQUAL( testcriterion3->getPrecision(), defaultEps );
    BOOST_CHECK_EQUAL( testcriterion3->getCheckMode(), ResidualCheck::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetCheckModeTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getCheckMode(), ResidualCheck::Absolute );
    mCriterionDouble->setCheckMode( ResidualCheck::Relative );
    BOOST_CHECK_EQUAL( mCriterionDouble->getCheckMode(), ResidualCheck::Relative );
    BOOST_CHECK_EQUAL( mCriterionFloat->getCheckMode(), ResidualCheck::Relative );
    mCriterionFloat->setCheckMode( ResidualCheck::Absolute );
    BOOST_CHECK_EQUAL( mCriterionFloat->getCheckMode(), ResidualCheck::Absolute );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetPrecisionTest )
{
    ValueType epsD  = 1e-6;
    ValueType epsD1 = 1e-5;
    BOOST_CHECK_EQUAL( mCriterionDouble->getPrecision(), epsD );
    mCriterionDouble->setPrecision( epsD1 );
    BOOST_CHECK_EQUAL( mCriterionDouble->getPrecision(), epsD1 );

    ValueType epsF  = 1e-4;
    ValueType epsF1 = 1e-3;
    BOOST_CHECK_EQUAL( mCriterionFloat->getPrecision(), epsF );
    mCriterionFloat->setPrecision( epsF1 );
    BOOST_CHECK_EQUAL( mCriterionFloat->getPrecision(), epsF1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    std::unique_ptr<ResidualThreshold<ValueType>> testcriterion1( mCriterionFloat->copy() );
    BOOST_CHECK_EQUAL( testcriterion1->getFirstNormResult(), -1.0 );
    BOOST_CHECK_EQUAL( testcriterion1->getPrecision(), ValueType( 1e-4) );
    BOOST_CHECK_EQUAL( testcriterion1->getCheckMode(), ResidualCheck::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetFirstNormResultTest )
{
    std::unique_ptr<ResidualThreshold<ValueType>> testcriterion( new ResidualThreshold<ValueType>() );
    BOOST_CHECK_EQUAL( testcriterion->getFirstNormResult(), ValueType(-1) );
    RealType<ValueType> s = 1;
    testcriterion->setFirstNormResult( s );
    BOOST_CHECK_EQUAL( testcriterion->getFirstNormResult(), s );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( NormPtrTest )
{
    NormPtr<ValueType> norm = mCriterionFloat->getNorm();

    // Check if this norm is really a L2Norm.

    DenseVector<ValueType> v( HArray<ValueType>( { 3, 4 } ) );

    RealType<ValueType> computed = norm->apply( v );
    RealType<ValueType> expected = 5;

    SCAI_LOG_INFO( logger, "v = " << v << ", Norm = " << *norm 
                           << ", computed = " << computed << ", expected = " << expected )

    BOOST_CHECK_CLOSE( computed, expected, 0.001 );
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
void testIsSatisfied( ResidualCheck checkMode )
{
    IndexType n = 3;
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();
    Jacobi<ValueType> solver( "ResidualThresholdTest Jacobi", ( ValueType ) 1.0 );
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );

    auto solution = denseVectorFill<ValueType>( n, 0 );

    auto l2Norm = std::make_shared<L2Norm<ValueType>>();

    RealType<ValueType> eps = 1e-5;
    CriterionPtr<ValueType> residualThreshold( new ResidualThreshold<ValueType>( l2Norm, eps, checkMode ) );
    solver.setStoppingCriterion( residualThreshold );
    solver.initialize( coefficients );
    solver.solve( solution, rhs );

    BOOST_CHECK( residualThreshold->isSatisfied( solver ) );

    const auto residual = denseVectorEval( coefficients * solution - rhs );

    RealType<ValueType> error = ( *l2Norm )( residual );

    SCAI_LOG_INFO( logger, "solver = " << solver << ", error = " << error )

    if ( checkMode == ResidualCheck::Absolute )
    {
        BOOST_CHECK( error <= eps );
    }
    else
    {
        BOOST_CHECK( error >= eps );        // relative comparison mode is less accurate
        BOOST_CHECK( error <= 10 * eps );   // but not this far
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, ValueType, scai_numeric_test_types )
{
    testIsSatisfied<ValueType>( ResidualCheck::Absolute );
    testIsSatisfied<ValueType>( ResidualCheck::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    std::unique_ptr<ResidualThreshold<ValueType>> testcriterion( new ResidualThreshold<ValueType>() );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
