/**
 * @file BlockDistributionTest.cpp
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
 * @brief Contains the implementation of the class BlockDistributionTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/distribution/BlockDistribution.hpp>

#include <scai/lama/test/distributed/DistributionTest.hpp>
#include <scai/lama/test/TestMacros.hpp>

using namespace scai::lama;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct BlockDistributionTestConfig
{
    BlockDistributionTestConfig()
    {
        comm = Communicator::getCommunicator( scai::lama::communicator::MPI );
        rank = comm->getRank();
        size = comm->getSize();
        blockSize = 17;
        dist = DistributionPtr( new BlockDistribution( blockSize * size, comm ) );
    }

    ~BlockDistributionTestConfig()
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

typedef boost::mpl::list<double, float> test_types;

BOOST_FIXTURE_TEST_SUITE( BlockDistributionTest, BlockDistributionTestConfig )
;

SCAI_LOG_DEF_LOGGER( logger, "Test.BlockDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    DistributionTest disttest( dist );

    if ( base_test_case )
    {
        SCAI_LOG_INFO( logger, "Run test method " << testcase << " in BlockDistributionTest." );
        DISTRIBUTION_COMMONTESTCASES( disttest );
    }
    else
    {
        disttest.runTests();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CreateTest )
{
    DistributionPtr bdist ( BlockDistribution::create( DistributionArguments( comm, 1, NULL, 1.0 ) ) );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), 1 );
    bdist.reset( Distribution::getDistribution( "BLOCK", comm, 1 ) );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), 1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( blockSizeTest )
{
    // Test smaller sizes
    for ( IndexType n = 1; n <= size; n++ )
    {
        BlockDistribution small( n, comm );

        // only the first n partitions have one element

        if ( rank < n )
        {
            BOOST_CHECK( small.getLocalSize() == 1 );
        }
        else
        {
            BOOST_CHECK( small.getLocalSize() == 0 );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( blockComputeOwnersTest )
{
    for ( PartitionId p = 0; p < size; ++p )
    {
        if ( p == rank )
        {
            for ( IndexType i = 0; i < blockSize; ++i )
            {
                BOOST_CHECK( dist->isLocal( p * blockSize + i ) );
                BOOST_CHECK_EQUAL( dist->global2local( p * blockSize + i ), i );
            }
        }
        else
        {
            for ( IndexType i = 0; i < blockSize; ++i )
            {
                nonLocalIndexes.push_back( p * blockSize + i );
                BOOST_CHECK_EQUAL( dist->global2local( p * blockSize + i ), nIndex );
            }
        }
    }

    std::vector<PartitionId> owners;
    dist->computeOwners( nonLocalIndexes, owners );
    BOOST_CHECK_EQUAL( ( int ) owners.size(), ( size - 1 ) * blockSize );
    std::vector<PartitionId>::size_type currentIndex = 0;

    for ( PartitionId p = 0; p < size; ++p )
    {
        if ( p == rank )
        {
            continue;
        }

        for ( IndexType i = 0; i < blockSize; ++i )
        {
            BOOST_CHECK_EQUAL( p, owners[currentIndex++] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    DistributionPtr blockdist1( new BlockDistribution( 1, comm ) );
    DistributionPtr blockdist2( blockdist1 );
    DistributionPtr blockdist3( new BlockDistribution( 1, comm ) );
    DistributionPtr blockdist4( new BlockDistribution( 2, comm ) );
    BOOST_CHECK( ( *blockdist1 ).isEqual( *blockdist2 ) );
    BOOST_CHECK( ( *blockdist1 ).isEqual( *blockdist3 ) );
    BOOST_CHECK( !( *blockdist1 ).isEqual( *blockdist4 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
