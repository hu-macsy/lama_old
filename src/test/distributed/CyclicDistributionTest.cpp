/**
 * @file CyclicDistributionTest.cpp
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
 * @brief Contains the implementation of the class CyclicDistributionTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <lama/distribution/CyclicDistribution.hpp>

#include <test/distributed/DistributionTest.hpp>
#include <test/TestMacros.hpp>

using namespace lama;
using namespace boost;

extern bool base_test_case;
extern std::string testcase;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct CyclicDistributionTestConfig
{
    CyclicDistributionTestConfig()
    {
        comm = Communicator::get( "MPI" );
        size = comm->getSize();
        rank = comm->getRank();
        chunkSize = size;
        globalSize = 10 * size;
        dist = DistributionPtr( new CyclicDistribution( globalSize, chunkSize, comm ) );
    }

    ~CyclicDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType chunkSize;
    IndexType globalSize;

    DistributionPtr dist;
};

BOOST_FIXTURE_TEST_SUITE( CyclicDistributionTest, CyclicDistributionTestConfig )

LAMA_LOG_DEF_LOGGER( logger, "Test.CyclicDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( commonTestCases )
{
    DistributionTest disttest( dist );

    if ( base_test_case )
    {
        LAMA_LOG_INFO( logger, "Run test method " << testcase << " in CyclicDistributionTest." );
        DISTRIBUTION_COMMONTESTCASES( disttest );
    }
    else
    {
        disttest.runTests();
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( cyclicGlobal2Local )
{
    //test local2global
    IndexType counter = rank * chunkSize;

    for ( IndexType i = 0; i < dist->getLocalSize(); i++ )
    {
        if ( i != 0 && i % chunkSize == 0 )
        {
            counter += ( size - 1 ) * chunkSize;
        }

        BOOST_CHECK_EQUAL( dist->local2global( i ), counter );
        counter++;
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( cyclicComputeOwnersTest )
{
    chunkSize = 3;
    globalSize = 16;
    std::vector<IndexType> localSizes;
    std::vector<PartitionId> theOwners; // used for verification
    CyclicDistribution distribution( globalSize, chunkSize, comm );
    IndexType chunkOffset = 0;
    IndexType chunkNr = 0;

    for ( IndexType i = 0; i < globalSize; i++ )
    {
        PartitionId owner = chunkNr % size;
        theOwners.push_back( owner );
        chunkOffset++;

        if ( chunkOffset == chunkSize )
        {
            chunkOffset = 0;
            chunkNr++;
        }
    }

    std::vector<IndexType> indexes;

    for ( IndexType i = 0; i < globalSize; i++ )
    {
        indexes.push_back( i );
    }

    std::vector<PartitionId> owners;
    distribution.computeOwners( indexes, owners );
    BOOST_CHECK_EQUAL( globalSize, static_cast<IndexType>( owners.size() ) );
    BOOST_CHECK_EQUAL( globalSize, static_cast<IndexType>( theOwners.size() ) );

    // now check for correct owners
    for ( IndexType i = 0; i < globalSize; i++ )
    {
        BOOST_CHECK_EQUAL( theOwners[i], owners[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( cyclicTest )
{
    // Test cyclic distribution for different global sizes and chunk sizes
    IndexType globalSizes[] =
    { 1, 50, 200 };
    IndexType chunkSizes[] =
    { 1, 3, 5, 19 };

    //BOOST_FOREACH(IndexType globalSize, globalSizes)
    for ( IndexType i = 0; i < 3; ++i )
    {
        IndexType globalSize = globalSizes[i];

        //BOOST_FOREACH(IndexType chunkSize, chunkSizes)
        for ( IndexType j = 0; j < 4; ++j )
        {
            IndexType chunkSize = chunkSizes[i];
            // printf("Test cyclic distribution of %d elements in chunks of size %d\n", globalSize, chunkSize);
            CyclicDistribution cyclicDist( globalSize, chunkSize, comm );
            // test the routines getNumChunks + getNumTotalChunks
            IndexType sumLocalChunks = cyclicDist.getCommunicator().sum( cyclicDist.getNumLocalChunks() );
            BOOST_CHECK_EQUAL( cyclicDist.getNumTotalChunks(), sumLocalChunks );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    DistributionPtr cyclicdist1( new CyclicDistribution( 1, comm->getSize(), comm ) );
    DistributionPtr cyclicdist2( cyclicdist1 );
    DistributionPtr cyclicdist3( new CyclicDistribution( 1, comm->getSize(), comm ) );
    DistributionPtr cyclicdist4( new CyclicDistribution( 3, comm->getSize(), comm ) );
    BOOST_CHECK( ( *cyclicdist1 ).isEqual( *cyclicdist2 ) );
    BOOST_CHECK( ( *cyclicdist1 ).isEqual( *cyclicdist3 ) );
    BOOST_CHECK( !( *cyclicdist1 ).isEqual( *cyclicdist4 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
