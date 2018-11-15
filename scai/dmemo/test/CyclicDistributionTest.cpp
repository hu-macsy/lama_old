/**
 * @file CyclicDistributionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2018
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
 *
 * LAMA is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * LAMA is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for
 * more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with LAMA. If not, see <http://www.gnu.org/licenses/>.
 * @endlicense
 *
 * @brief Contains the implementation of the class CyclicDistributionTest.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/utilskernel/HArrayUtils.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct CyclicDistributionTestConfig
{
    CyclicDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
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

SCAI_LOG_DEF_LOGGER( logger, "Test.CyclicDistributionTest" );

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
        theOwners.push_back( chunkNr );

        chunkOffset++;

        if ( chunkOffset == chunkSize )
        {
            chunkOffset = 0;
            chunkNr++;

            if ( chunkNr == static_cast<IndexType>( size ) )
            {
                chunkNr = 0;
            }
        }
    }

    using namespace utilskernel;

    hmemo::HArray<IndexType> indexes;
    utilskernel::HArrayUtils::setOrder( indexes, globalSize );

    hmemo::HArray<PartitionId> owners;
    distribution.computeOwners( owners, indexes );

    BOOST_CHECK_EQUAL( owners.size(), indexes.size() );

    hmemo::ReadAccess<PartitionId> rOwners( owners );

    BOOST_CHECK_EQUAL( globalSize, owners.size() );
    BOOST_CHECK_EQUAL( globalSize, static_cast<IndexType>( theOwners.size() ) );

    // now check for correct owners

    for ( IndexType i = 0; i < globalSize; i++ )
    {
        BOOST_CHECK_EQUAL( theOwners[i], rOwners[i] );
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
            IndexType sumLocalChunks = cyclicDist.getTargetCommunicator().sum( cyclicDist.getNumLocalChunks() );
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
