/**
 * @file GenBlockDistributionTest.cpp
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
 * @brief Contains the implementation of the class GenBlockDistributionTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/lama/distribution/GenBlockDistribution.hpp>

#include <test/distributed/DistributionTest.hpp>
#include <test/TestMacros.hpp>

using namespace lama;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct GenBlockDistributionTestConfig
{
    GenBlockDistributionTestConfig()
    {
        comm = Communicator::get( "MPI" );
        size = comm->getSize();
        rank = comm->getRank();
        globalSize = size * ( size + 1 );

        for ( PartitionId p = 0; p < size; ++p )
        {
            IndexType localSize = 2 * ( p + 1 );
            localSizes.push_back( localSize );

            for ( IndexType i = 0; i < localSize; i++ )
            {
                theOwners.push_back( p );
            }
        }

        dist = DistributionPtr( new GenBlockDistribution( globalSize, localSizes, comm ) );
    }

    ~GenBlockDistributionTestConfig()
    {
        comm = CommunicatorPtr();
        dist = DistributionPtr();
    }

    DistributionPtr dist;

    PartitionId size;
    PartitionId rank;
    IndexType globalSize;

    std::vector<IndexType> localSizes;
    std::vector<PartitionId> theOwners;
};

BOOST_FIXTURE_TEST_SUITE( GenBlockDistributionTest, GenBlockDistributionTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.GenBlockDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    DistributionTest disttest( dist );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run test method " << testcase << " in GenBlockDistributionTest." );
        DISTRIBUTION_COMMONTESTCASES( disttest );
    }
    else
    {
        disttest.runTests();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( genBlockComputeOwnersTest )
{
    std::vector<IndexType> indexes;

    for ( IndexType i = 0; i < globalSize; i++ )
    {
        indexes.push_back( i );
    }

    std::vector<PartitionId> owners;
    dist->computeOwners( indexes, owners );
    BOOST_CHECK_EQUAL( globalSize, static_cast<IndexType>( owners.size() ) );
    BOOST_CHECK_EQUAL( globalSize, static_cast<IndexType>( theOwners.size() ) );

    // now check for correct owners

    for ( IndexType i = 0; i < globalSize; i++ )
    {
        BOOST_CHECK_EQUAL( theOwners[i], owners[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( genBlockSizeTest )
{
    // setup vector of sizes
    std::vector<IndexType> mlocalSizes; // block size for each partition

    for ( PartitionId p = 0; p < size; ++p )
    {
        // Partition p gets '(p + 1) * 2' elements
        mlocalSizes.push_back( 2 * ( p + 1 ) );
    }

    IndexType globalSize = size * ( size + 1 );
    GenBlockDistribution dist1( globalSize, mlocalSizes, comm );
    BOOST_CHECK( dist1.getLocalSize() == 2 * ( rank + 1 ) );
    IndexType lb1, ub1;
    dist1.getLocalRange( lb1, ub1 );
    GenBlockDistribution dist2( globalSize, 2 * ( rank + 1 ), comm );
    IndexType lb2, ub2;
    dist1.getLocalRange( lb2, ub2 );
    BOOST_CHECK( lb1 == lb2 );
    BOOST_CHECK( ub1 == ub2 );
    GenBlockDistribution dist3( globalSize, lb1, ub1, comm );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    DistributionPtr genblockdist1( new GenBlockDistribution( comm->getSize(), 1, comm ) );
    DistributionPtr genblockdist2( genblockdist1 );
    DistributionPtr genblockdist3( new GenBlockDistribution( comm->getSize(), 1, comm ) );
    DistributionPtr genblockdist4( new GenBlockDistribution( 2 * comm->getSize(), 2, comm ) );
    BOOST_CHECK( ( *genblockdist1 ).isEqual( *genblockdist2 ) );
    BOOST_CHECK( ( *genblockdist1 ).isEqual( *genblockdist3 ) );
    BOOST_CHECK( !( *genblockdist1 ).isEqual( *genblockdist4 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
