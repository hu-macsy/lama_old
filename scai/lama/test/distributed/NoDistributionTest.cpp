/**
 * @file NoDistributionTest.cpp
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
 * @brief Contains the implementation of the class NoDistributionTest.
 * @author: Alexander Büchel
 * @date 01.08.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/distribution/NoDistribution.hpp>

#include <test/distributed/DistributionTest.hpp>
#include <test/TestMacros.hpp>

using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct NoDistributionTestConfig
{
    NoDistributionTestConfig()
    {
        comm = Communicator::get( "MPI" );
        rank = comm->getRank();
        size = comm->getSize();
        blockSize = 17;
        dist = DistributionPtr( new NoDistribution( blockSize * size ) );
    }

    ~NoDistributionTestConfig()
    {
        comm = CommunicatorPtr();
        dist = DistributionPtr();
    }

    PartitionId rank;
    PartitionId size;
    IndexType blockSize;

    std::vector<IndexType> nonLocalIndexes;

    DistributionPtr dist;
};

BOOST_FIXTURE_TEST_SUITE( NoDistributionTest, NoDistributionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.NoDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    DistributionPtr nodist( new NoDistribution( 20 ) );
    BOOST_CHECK_EQUAL( nodist->getCommunicatorPtr()->getType(), "none" );
    BOOST_CHECK_EQUAL( nodist->getGlobalSize(), 20 );
    BOOST_CHECK_EQUAL( nodist->getLocalSize(), 20 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    DistributionTest disttest( dist );

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run test method " << testcase << " in NoDistributionTest." );
        DISTRIBUTION_COMMONTESTCASES( disttest );
    }
    else
    {
        disttest.runTests();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isLocalTest )
{
    BOOST_CHECK( dist->isLocal( 0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    DistributionPtr nodist1( new NoDistribution( 1 ) );
    DistributionPtr nodist2( nodist1 );
    DistributionPtr nodist3( new NoDistribution( 1 ) );
    DistributionPtr nodist4( new NoDistribution( 2 ) );
    BOOST_CHECK( ( *nodist1 ).isEqual( *nodist2 ) );
    BOOST_CHECK( ( *nodist1 ).isEqual( *nodist3 ) );
    BOOST_CHECK( !( *nodist1 ).isEqual( *nodist4 ) );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
