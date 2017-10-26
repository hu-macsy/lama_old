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
#include <scai/common/unique_ptr.hpp>
#include <scai/solver/CG.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai;
using namespace solver;
using namespace hmemo;

/* --------------------------------------------------------------------- */

struct CriterionTestConfig
{
    CriterionTestConfig()
    {
        mIterationCountCriterion1 = new IterationCount( 10 );
        mIterationCountCriterion1Ptr = CriterionPtr( mIterationCountCriterion1 );
        mIterationCountCriterion2 = new IterationCount( 8 );
        mIterationCountCriterion2Ptr = CriterionPtr( mIterationCountCriterion2 );
    }

    ~CriterionTestConfig()
    {
    }

    CriterionPtr mCriterionRoot;
    IterationCount* mIterationCountCriterion1;
    CriterionPtr mIterationCountCriterion1Ptr;
    IterationCount* mIterationCountCriterion2;
    CriterionPtr mIterationCountCriterion2Ptr;
};

BOOST_FIXTURE_TEST_SUITE( CriterionTest, CriterionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.CriterionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mIterationCountCriterion1->getIterationExtrema(), IndexType( 10 ) );
    BOOST_CHECK_EQUAL( mIterationCountCriterion2->getIterationExtrema(), IndexType( 8 ) );
    Criterion* boolcondition1 = new Criterion();
    CriterionPtr testcriterion1 = CriterionPtr( boolcondition1 );
    BOOST_CHECK_EQUAL( testcriterion1->getOperation(), Criterion::AND );
    Criterion* boolcondition2 = new Criterion( true );
    CriterionPtr testcriterion2 = CriterionPtr( boolcondition2 );
    BOOST_CHECK_EQUAL( testcriterion2->getOperation(), Criterion::AND );
    Criterion* boolcondition3 = new Criterion( boolcondition1 );
    CriterionPtr testcriterion3 = CriterionPtr( boolcondition3 );
    BOOST_CHECK_EQUAL( testcriterion3->getOperation(), testcriterion1->getOperation() );
    BOOST_CHECK_EQUAL( testcriterion3->getLeftChild(), testcriterion1->getLeftChild() );
    BOOST_CHECK_EQUAL( testcriterion3->getRightChild(), testcriterion1->getRightChild() );
    Criterion* boolcondition4 = new Criterion( testcriterion1, true );
    CriterionPtr testcriterion4 = CriterionPtr( boolcondition4 );
    BOOST_CHECK_EQUAL( testcriterion4->getOperation(), testcriterion1->getOperation() );
    BOOST_CHECK_EQUAL( testcriterion4->getLeftChild(), testcriterion1->getLeftChild() );
    BOOST_CHECK_EQUAL( testcriterion4->getRightChild(), testcriterion1->getRightChild() );
    Criterion* boolcondition5 = new Criterion( mIterationCountCriterion1Ptr, mIterationCountCriterion2Ptr,
            Criterion::AND );
    CriterionPtr testcriterion5 = CriterionPtr( boolcondition5 );
    BOOST_CHECK_EQUAL( testcriterion5->getOperation(), Criterion::AND );
    BOOST_CHECK_EQUAL( testcriterion5->getLeftChild(), mIterationCountCriterion1Ptr );
    BOOST_CHECK_EQUAL( testcriterion5->getRightChild(), mIterationCountCriterion2Ptr );
}

/* --------------------------------------------------------------------- */

void testCheckCriterion( CriterionPtr bcp, Criterion::BooleanOperator boperator )
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
    mCriterionRoot = mIterationCountCriterion1Ptr && mIterationCountCriterion2Ptr;
    testCheckCriterion( mCriterionRoot, Criterion::AND );
    mCriterionRoot = mIterationCountCriterion1Ptr || mIterationCountCriterion2Ptr;
    testCheckCriterion( mCriterionRoot, Criterion::OR );
    //TODO negation!!!
    //mCriterionRoot = !(mIterationCountCriterion1Ptr || mIterationCountCriterion2Ptr) && ((mIterationCountCriterion2Ptr && !true) || false);
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
    Criterion* boolcondition = new Criterion();
    CriterionPtr testcriterion = CriterionPtr( boolcondition );
    testcriterion->setLeftChild( mIterationCountCriterion1Ptr );
    testcriterion->setRightChild( mIterationCountCriterion2Ptr );
    testcriterion->setOperation( Criterion::OR );
    BOOST_CHECK_EQUAL( testcriterion->getLeftChild(), mIterationCountCriterion1Ptr );
    BOOST_CHECK_EQUAL( testcriterion->getRightChild(), mIterationCountCriterion2Ptr );
    BOOST_CHECK_EQUAL( testcriterion->getOperation(), Criterion::OR );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    Criterion* boolcondition = new Criterion();
    CriterionPtr testcriterion = CriterionPtr( boolcondition );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion ); // no empty string for mModifier = true
    Criterion* boolcondition2 = new Criterion( false );
    CriterionPtr testcriterion2 = CriterionPtr( boolcondition2 );
    SCAI_COMMON_WRITEAT_TEST( *testcriterion2 ); // no empty string for mModifier = false
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE ( isSatisfiedTest )
{
    typedef SCAI_TEST_TYPE ValueType;
    const IndexType N = 40;
    scai::lama::CSRSparseMatrix<ValueType> coefficients;
    scai::lama::MatrixCreator::buildPoisson2D( coefficients, 5, N, N );
    scai::lama::DenseVector<ValueType> rhs( coefficients.getRowDistributionPtr(), 1.0 );
    scai::lama::DenseVector<ValueType> solution( coefficients.getColDistributionPtr(), 1.0 );
    CG cgsolver( "CriterionTestSolver" );
    cgsolver.setStoppingCriterion( mIterationCountCriterion2Ptr );
    cgsolver.initialize( coefficients );
    cgsolver.solve( solution, rhs );
    bool test;
    // isSatisfied for BooleanCondition ( noChildren, returns mModifier )
    common::unique_ptr<Criterion> condition1( new Criterion() );
    test = condition1->isSatisfied( cgsolver );
    BOOST_CHECK ( test );
    // isSatisfied for BooleanCondition ( noChildren, returns mModifier )
    common::unique_ptr<Criterion> condition2( new Criterion( false ) );
    test = condition2->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left
    common::unique_ptr<Criterion> condition3( new Criterion( false ) );
    condition3->setLeftChild( mIterationCountCriterion1Ptr );
    test = condition3->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left
    common::unique_ptr<Criterion> condition4( new Criterion( false ) );
    condition4->setRightChild( mIterationCountCriterion1Ptr );
    test = condition4->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left AND right child
    common::unique_ptr<Criterion> condition5( new Criterion( mIterationCountCriterion1Ptr, mIterationCountCriterion2Ptr, Criterion::AND ) );
    test = condition5->isSatisfied( cgsolver );
    BOOST_CHECK ( !test );
    // isSatisfied for Criterion with left OR right child
    common::unique_ptr<Criterion> condition6( new Criterion( mIterationCountCriterion1Ptr, mIterationCountCriterion2Ptr, Criterion::OR ) );
    test = condition6->isSatisfied( cgsolver );
    BOOST_CHECK ( test );
}

/* ---------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
