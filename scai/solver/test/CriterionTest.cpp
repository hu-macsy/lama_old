/**
 * @file CriterionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 * @brief CriterionTest.hpp
 * @author Kai Buschulte
 * @date 02.02.2012
 */

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <scai/solver/criteria/Criterion.hpp>
#include <scai/solver/criteria/IterationCount.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/matrix/CSRSparseMatrix.hpp>
#include <scai/lama/matutils/MatrixCreator.hpp>
#include <scai/solver/CG.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai;
using namespace scai::solver;
using namespace scai::hmemo;

typedef SCAI_TEST_TYPE ValueType;

/* --------------------------------------------------------------------- */

struct CriterionTestConfig
{
    typedef SCAI_TEST_TYPE ValueType;

    CriterionTestConfig()
    {
        mIterationCountCriterion1 = new IterationCount<ValueType>( 10 );
        mIterationCountCriterion1Ptr.reset( mIterationCountCriterion1 );
        mIterationCountCriterion2 = new IterationCount<ValueType>( 8 );
        mIterationCountCriterion2Ptr.reset( mIterationCountCriterion2 );
    }

    ~CriterionTestConfig()
    {
    }

    IterationCount<ValueType>* mIterationCountCriterion1;
    CriterionPtr<ValueType> mIterationCountCriterion1Ptr;
    IterationCount<ValueType>* mIterationCountCriterion2;
    CriterionPtr<ValueType> mIterationCountCriterion2Ptr;
};

BOOST_FIXTURE_TEST_SUITE( CriterionTest, CriterionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.CriterionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mIterationCountCriterion1->getIterationExtrema(), IndexType( 10 ) );
    BOOST_CHECK_EQUAL( mIterationCountCriterion2->getIterationExtrema(), IndexType( 8 ) );
    Criterion<ValueType>* boolcondition1 = new Criterion<ValueType>();
    CriterionPtr<ValueType> testcriterion1( boolcondition1 );
    BOOST_CHECK_EQUAL( testcriterion1->getOperation(), BooleanOp::AND );
    Criterion<ValueType>* boolcondition2 = new Criterion<ValueType>( true );
    CriterionPtr<ValueType> testcriterion2 = CriterionPtr<ValueType>( boolcondition2 );
    BOOST_CHECK_EQUAL( testcriterion2->getOperation(), BooleanOp::AND );
    Criterion<ValueType>* boolcondition3 = new Criterion<ValueType>( boolcondition1 );
    CriterionPtr<ValueType> testcriterion3 = CriterionPtr<ValueType>( boolcondition3 );
    BOOST_CHECK_EQUAL( testcriterion3->getOperation(), testcriterion1->getOperation() );
    BOOST_CHECK_EQUAL( testcriterion3->getLeftChild(), testcriterion1->getLeftChild() );
    BOOST_CHECK_EQUAL( testcriterion3->getRightChild(), testcriterion1->getRightChild() );
    Criterion<ValueType>* boolcondition4 = new Criterion<ValueType>( testcriterion1, true );
    CriterionPtr<ValueType> testcriterion4 = CriterionPtr<ValueType>( boolcondition4 );
    BOOST_CHECK_EQUAL( testcriterion4->getOperation(), testcriterion1->getOperation() );
    BOOST_CHECK_EQUAL( testcriterion4->getLeftChild(), testcriterion1->getLeftChild() );
    BOOST_CHECK_EQUAL( testcriterion4->getRightChild(), testcriterion1->getRightChild() );
    Criterion<ValueType>* boolcondition5 = new Criterion<ValueType>( mIterationCountCriterion1Ptr, mIterationCountCriterion2Ptr,
            BooleanOp::AND );
    CriterionPtr<ValueType> testcriterion5( boolcondition5 );
    BOOST_CHECK_EQUAL( testcriterion5->getOperation(), BooleanOp::AND );
    BOOST_CHECK_EQUAL( testcriterion5->getLeftChild(), mIterationCountCriterion1Ptr );
    BOOST_CHECK_EQUAL( testcriterion5->getRightChild(), mIterationCountCriterion2Ptr );
}

/* --------------------------------------------------------------------- */

template<typename ValueType>
void testCheckCriterion( CriterionPtr<ValueType> bcp, BooleanOp boperator )
{
    BOOST_CHECK_EQUAL( bcp->getOperation(), boperator );
    BOOST_CHECK( bcp->hasLeftChild() );
    BOOST_CHECK( bcp->hasRightChild() );
    BOOST_CHECK( !bcp->getLeftChild()->hasLeftChild() );
    BOOST_CHECK( !bcp->getRightChild()->hasLeftChild() );
    BOOST_CHECK( !bcp->getLeftChild()->hasRightChild() );
    BOOST_CHECK( !bcp->getRightChild()->hasRightChild() );
}

BOOST_AUTO_TEST_CASE( ConditionTest )
{
    CriterionPtr<ValueType> mCriterionRoot;
    mCriterionRoot = mIterationCountCriterion1Ptr && mIterationCountCriterion2Ptr;
    testCheckCriterion( mCriterionRoot, BooleanOp::AND );
    mCriterionRoot = mIterationCountCriterion1Ptr || mIterationCountCriterion2Ptr;
    testCheckCriterion( mCriterionRoot, BooleanOp::OR );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( AssignmentTest )
{
    mIterationCountCriterion1 = mIterationCountCriterion2;
    BOOST_CHECK_EQUAL( mIterationCountCriterion1->getIterationExtrema(), IndexType( 8 ) );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( SetAndGetChildTest )
{
    Criterion<ValueType>* boolcondition = new Criterion<ValueType>();
    CriterionPtr<ValueType> testcriterion = CriterionPtr<ValueType>( boolcondition );
    testcriterion->setLeftChild( mIterationCountCriterion1Ptr );
    testcriterion->setRightChild( mIterationCountCriterion2Ptr );
    testcriterion->setOperation( BooleanOp::OR );
    BOOST_CHECK_EQUAL( testcriterion->getLeftChild(), mIterationCountCriterion1Ptr );
    BOOST_CHECK_EQUAL( testcriterion->getRightChild(), mIterationCountCriterion2Ptr );
    BOOST_CHECK_EQUAL( testcriterion->getOperation(), BooleanOp::OR );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    Criterion<ValueType>* boolcondition = new Criterion<ValueType>();
    CriterionPtr<ValueType> testcriterion = CriterionPtr<ValueType>( boolcondition );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion ); // no empty string for mModifier = true
    Criterion<ValueType>* boolcondition2 = new Criterion<ValueType>( false );
    CriterionPtr<ValueType> testcriterion2 = CriterionPtr<ValueType>( boolcondition2 );
    // SCAI_COMMON_WRITEAT_TEST( *testcriterion2 ); // no empty string for mModifier = false
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( isSatisfiedTest )
{
    const IndexType N = 40;

    lama::CSRSparseMatrix<ValueType> coefficients;
    lama::MatrixCreator::buildPoisson2D( coefficients, 5, N, N );

    auto rhs      = lama::fill<lama::DenseVector<ValueType>>( coefficients.getRowDistributionPtr(), ValueType( 1 ) );
    auto solution = lama::fill<lama::DenseVector<ValueType>>( coefficients.getColDistributionPtr(), ValueType( 1 ) );

    CG<ValueType> cgsolver( "CriterionTestSolver" );
    cgsolver.setStoppingCriterion( mIterationCountCriterion2Ptr );
    cgsolver.initialize( coefficients );
    cgsolver.solve( solution, rhs );
    bool test;
    // isSatisfied for BooleanCondition ( noChildren, returns mModifier )
    Criterion<ValueType>* condition1 = new Criterion<ValueType>();
    test = condition1->isSatisfied( cgsolver );
    BOOST_CHECK ( test );
    // isSatisfied for BooleanCondition ( noChildren, returns mModifier )
    Criterion<ValueType>* condition2 = new Criterion<ValueType>( false );
    test = condition2->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left
    Criterion<ValueType>* condition3 = new Criterion<ValueType>( false );
    condition3->setLeftChild( mIterationCountCriterion1Ptr );
    test = condition3->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left
    Criterion<ValueType>* condition4 = new Criterion<ValueType>( false );
    condition4->setRightChild( mIterationCountCriterion1Ptr );
    test = condition4->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left AND right child
    Criterion<ValueType>* condition5 = new Criterion<ValueType>( mIterationCountCriterion1Ptr, mIterationCountCriterion2Ptr, BooleanOp::AND );
    test = condition5->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left OR right child
    Criterion<ValueType>* condition6 = new Criterion<ValueType>( mIterationCountCriterion1Ptr, mIterationCountCriterion2Ptr, BooleanOp::OR );
    test = condition6->isSatisfied( cgsolver );
    BOOST_CHECK ( test );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
