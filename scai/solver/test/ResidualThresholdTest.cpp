/**
 * @file ResidualThresholdTest.cpp
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
 * @brief Contains the implementation of the class ResidualThresholdTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 14.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/Scalar.hpp>
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

#include <scai/common/unique_ptr.hpp>

using namespace scai;
using namespace scai::solver;
using namespace scai::lama;
using namespace scai::hmemo;

/* --------------------------------------------------------------------- */

struct ResidualThresholdTestConfig
{
    ResidualThresholdTestConfig()
    {
        NormPtr mNorm = NormPtr( new MaxNorm() );
        mCriterionDouble = new ResidualThreshold( mNorm, Scalar( 1.0e-6 ), ResidualThreshold::Absolute );
        mCriterionFloat = new ResidualThreshold( mNorm, Scalar( 1.0e-4f ), ResidualThreshold::Relative );
    }

    ~ResidualThresholdTestConfig()
    {
    }

    ResidualThreshold* mCriterionDouble;
    ResidualThreshold* mCriterionFloat;
};

BOOST_FIXTURE_TEST_SUITE( ResidualThresholdTest, ResidualThresholdTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.ResidualThresholdTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getPrecision().getValue<double>(), 1.0e-6 );
    BOOST_CHECK_EQUAL( mCriterionDouble->getCheckMode(), ResidualThreshold::Absolute );
    BOOST_CHECK_EQUAL( mCriterionFloat->getPrecision().getValue<float>(), 1.0e-4f );
    BOOST_CHECK_EQUAL( mCriterionFloat->getCheckMode(), ResidualThreshold::Relative );
    const ResidualThreshold* testcriterion1 = new ResidualThreshold();
    BOOST_CHECK_EQUAL( testcriterion1->getPrecision().getValue<float>(), 1.0e-5f );
    BOOST_CHECK_EQUAL( testcriterion1->getCheckMode(), ResidualThreshold::Relative );
    BOOST_CHECK_EQUAL( testcriterion1->getFirstNormResult(), -1.0 );
    NormPtr testnorm = NormPtr( new MaxNorm() );
    ResidualThreshold* testcriterion2 = new ResidualThreshold( testnorm );
    BOOST_CHECK_EQUAL( testcriterion2->getFirstNormResult(), -1.0 );
    BOOST_CHECK_EQUAL( testcriterion2->getPrecision().getValue<float>(), 1.0e-5f );
    BOOST_CHECK_EQUAL( testcriterion2->getCheckMode(), ResidualThreshold::Relative );
    ResidualThreshold* testcriterion3 = new ResidualThreshold( *testcriterion1 );
    BOOST_CHECK_EQUAL( testcriterion3->getFirstNormResult(), -1.0 );
    BOOST_CHECK_CLOSE( testcriterion3->getPrecision().getValue<float>(), 1.0e-5f, 1 );
    BOOST_CHECK_EQUAL( testcriterion3->getCheckMode(), ResidualThreshold::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetCheckModeTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getCheckMode(), ResidualThreshold::Absolute );
    mCriterionDouble->setCheckMode( ResidualThreshold::Relative );
    BOOST_CHECK_EQUAL( mCriterionDouble->getCheckMode(), ResidualThreshold::Relative );
    BOOST_CHECK_EQUAL( mCriterionFloat->getCheckMode(), ResidualThreshold::Relative );
    mCriterionFloat->setCheckMode( ResidualThreshold::Absolute );
    BOOST_CHECK_EQUAL( mCriterionFloat->getCheckMode(), ResidualThreshold::Absolute );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetPrecisionTest )
{
    BOOST_CHECK_EQUAL( mCriterionDouble->getPrecision().getValue<double>(), 1.0e-6 );
    mCriterionDouble->setPrecision( 1.0e-5 );
    BOOST_CHECK_EQUAL( mCriterionDouble->getPrecision().getValue<double>(), 1.0e-5 );
    BOOST_CHECK_EQUAL( mCriterionFloat->getPrecision().getValue<float>(), 1.0e-4f );
    mCriterionFloat->setPrecision( 1.0e-3f );
    BOOST_CHECK_EQUAL( mCriterionFloat->getPrecision().getValue<float>(), 1.0e-3f );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    ResidualThreshold* testcriterion1;
    testcriterion1 = ( ResidualThreshold* ) mCriterionFloat->copy();
    BOOST_CHECK_EQUAL( testcriterion1->getFirstNormResult(), -1.0 );
    BOOST_CHECK_CLOSE( testcriterion1->getPrecision().getValue<float>(), 1.0e-4f, 2 );
    BOOST_CHECK_EQUAL( testcriterion1->getCheckMode(), ResidualThreshold::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetFirstNormResultTest )
{
    Scalar s = 1.0;
    ResidualThreshold* testcriterion = new ResidualThreshold();
    BOOST_CHECK_EQUAL( testcriterion->getFirstNormResult(), -1.0 );
    testcriterion->setFirstNormResult( s );
    BOOST_CHECK_EQUAL( testcriterion->getFirstNormResult(), 1.0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( NormPtrTest )
{
    NormPtr norm = mCriterionFloat->getNorm();
    //TODO: Check if this norm is really a L2Norm.
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
void testIsSatisfied( ResidualThreshold::ResidualThresholdCheckMode checkMode )
{
    IndexType n = 3;
    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();
    Jacobi solver( "ResidualThresholdTest Jacobi", ( ValueType ) 1.0 );
    const CSRSparseMatrix<ValueType> coefficients( system.coefficients );
    const DenseVector<ValueType> rhs( system.rhs );
    DenseVector<ValueType> solution( n, 0.0 );
    NormPtr l2Norm( new L2Norm() );
    CriterionPtr residualThreshold( new ResidualThreshold( l2Norm, Scalar( 1e-6 ), checkMode ) );
    solver.setStoppingCriterion( residualThreshold );
    solver.initialize( coefficients );
    solver.solve( solution, rhs );
    DenseVector<ValueType> temp( n, 1.0 );
    temp = coefficients * solution;
    DenseVector<ValueType> error( n, 1.0 );
    error = rhs - temp;
    Scalar l2 = ( *l2Norm )( error );
    typedef typename common::TypeTraits<ValueType>::AbsType AbsType;
    BOOST_CHECK( l2.getValue<AbsType>() <= 1e-5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, ValueType, scai_arithmetic_test_types )
{
    testIsSatisfied<ValueType>( ResidualThreshold::Absolute );
    testIsSatisfied<ValueType>( ResidualThreshold::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    common::unique_ptr<ResidualThreshold> testcriterion( new ResidualThreshold() );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
