/**
 * @file P_NormTest.cpp
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
 * @brief Test norm for distributed vectors.
 * @author: Thomas Brandes
 * @date 10.05.2013
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/DenseVector.hpp>
#include <scai/lama/distribution/BlockDistribution.hpp>

#include <scai/lama/test/EquationHelper.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;
using namespace scai::hmemo;

typedef boost::mpl::list<float, double> test_types;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct P_NormTestConfig
{
    P_NormTestConfig()
    {
        comm = Communicator::getCommunicator();
    }

    ~P_NormTestConfig()
    {
        comm = CommunicatorPtr();
    }
};

/* --------------------------------------------------------------------- */

BOOST_FIXTURE_TEST_SUITE( P_NormTest, P_NormTestConfig );

SCAI_LOG_DEF_LOGGER( logger, "Test.P_NormTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( Norm, ValueType, test_types )
{
    CONTEXTLOOP()
    {
        GETCONTEXT( context );

        for ( IndexType size = 0; size < 4; ++size )
        {
            DistributionPtr dist( new BlockDistribution( size, comm ) );
            const ValueType VAL = 1.0;
            DenseVector<ValueType> repVector( size, VAL );
            DenseVector<ValueType> distVector;
            distVector.setContextPtr( context );
            distVector = repVector;
            distVector.redistribute( dist );
            Scalar l1norm = distVector.l1Norm();
            ValueType expectedL1Norm = static_cast<ValueType>( size );
            BOOST_CHECK_CLOSE( l1norm.getValue<ValueType>(), expectedL1Norm, 1 );
            Scalar l2norm = distVector.l2Norm();
            ValueType expectedL2Norm = static_cast<ValueType>( size );
            expectedL2Norm = std::sqrt( expectedL2Norm );
            BOOST_CHECK_CLOSE( l2norm.getValue<ValueType>(), expectedL2Norm, 1 );
            Scalar maxNorm = distVector.maxNorm();
            ValueType expectedMaxNorm = static_cast<ValueType>( VAL );

            if ( size == 0 )
            {
                expectedMaxNorm = static_cast<ValueType>( 0 );
            }

            BOOST_CHECK_CLOSE( maxNorm.getValue<ValueType>(), expectedMaxNorm, 1 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
