/**
 * @file ResidualThresholdTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class ResidualThresholdTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 14.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/Scalar.hpp>
#include <lama/DenseVector.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>

#include <lama/norm/MaxNorm.hpp>

#include <lama/solver/criteria/ResidualThreshold.hpp>
#include <lama/solver/DefaultJacobi.hpp>
#include <lama/norm/L2Norm.hpp>

#include <lama/expression/MatrixVectorExpressions.hpp>
#include <lama/expression/VectorExpressions.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace lama;
using namespace boost;

typedef boost::mpl::list<double,float> test_types;

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

LAMA_LOG_DEF_LOGGER( logger, "Test.ResidualThresholdTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ContstructorTest )
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
    testcriterion1 = (ResidualThreshold*) mCriterionFloat->copy();

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

template<typename DataType>
void testIsSatisfied( ResidualThreshold::ResidualThresholdCheckMode checkMode )
{
    typedef DataType ValueType;

    IndexType n = 3;

    EquationHelper::EquationSystem<ValueType> system = EquationHelper::get3x3SystemA<ValueType>();

    DefaultJacobi solver( "ResidualThresholdTest Jacobi", (ValueType) 1.0 );

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

    BOOST_CHECK( l2.getValue<ValueType>() <= 1e-5 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( redistributeTest, T, test_types ) {
    typedef T ValueType;

    testIsSatisfied<ValueType>( ResidualThreshold::Absolute );
    testIsSatisfied<ValueType>( ResidualThreshold::Relative );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    ResidualThreshold* testcriterion = new ResidualThreshold();

    LAMA_WRITEAT_PTR_TEST( testcriterion );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
