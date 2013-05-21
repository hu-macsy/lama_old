/**
 * @file IterationCountTest.cpp
 *
 * @license
 * Copyright (c) 2011
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
 * @brief Contains the implementation of the class IterationCountTest.
 * @author: Alexander Büchel, Matthias Makulla
 * @date 14.02.2012
 * $Id$
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/CommunicatorFactory.hpp>

#include <lama/solver/criteria/IterationCount.hpp>
#include <lama/solver/DefaultJacobi.hpp>
#include <lama/DenseVector.hpp>

#include <lama/matrix/CSRSparseMatrix.hpp>

#include <test/EquationHelper.hpp>
#include <test/TestMacros.hpp>

using namespace lama;
using namespace boost;

typedef boost::mpl::list<double,float> test_types;

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
;

LAMA_LOG_DEF_LOGGER( logger, "Test.IterationCountTest" );

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

template<typename DataType>
void testIsSatisfiedmethod( IndexType expectedIterationCount )
{
    typedef DataType ValueType;

    DefaultJacobi jacobi( "IterationCountTest Jacobi" );

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

BOOST_AUTO_TEST_CASE_TEMPLATE( testIsSatisfied, T, test_types ) {
    typedef T ValueType;

    testIsSatisfiedmethod<ValueType>( 10 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    LAMA_WRITEAT_PTR_TEST( mCriterionDouble );
    LAMA_WRITEAT_PTR_TEST( mCriterionFloat );

    std::stringstream mStream;
    mStream << mCriterionDouble;
    std::string mString = mStream.str();
    BOOST_CHECK( mString.compare("Maximal") );

    mStream.clear();
    mStream << mCriterionFloat;
    mString = mStream.str();
    BOOST_CHECK( mString.compare("Minimal") );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
