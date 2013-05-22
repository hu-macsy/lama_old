/**
 * @file CriterionTest.hpp
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
 * @brief CriterionTest.hpp
 * @author Alexander Büchel, Kai Buschulte
 * @date 02.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>
#include <iostream>

#include <lama/solver/criteria/Criterion.hpp>
#include <lama/solver/criteria/IterationCount.hpp>

#include <test/TestMacros.hpp>

using namespace boost;
using namespace lama;

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
;

LAMA_LOG_DEF_LOGGER( logger, "Test.CriterionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    BOOST_CHECK_EQUAL( mIterationCountCriterion1->getIterationExtrema(), 10 );

    BOOST_CHECK_EQUAL( mIterationCountCriterion2->getIterationExtrema(), 8 );

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
    BOOST_CHECK_EQUAL( mIterationCountCriterion1->getIterationExtrema(), 8 );
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
    LAMA_WRITEAT_PTR_TEST( testcriterion );
}
/* ---------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
