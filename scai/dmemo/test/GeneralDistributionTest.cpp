/**
 * @file GeneralDistributionTest.cpp
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
 * @brief Specific tests for derived distribution class GeneralDistributionTest.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/utilskernel.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct GeneralDistributionTestConfig
{
    GeneralDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
        rank = comm->getRank();
        size = comm->getSize();
        elemsPerPartition = 10;
        globalSize = elemsPerPartition * size;

        // get: rank, rank + size, rank + 2 * size, ...

        const IndexType first = static_cast<IndexType>( rank );
        const IndexType inc   = static_cast<IndexType>( size );

        utilskernel::HArrayUtils::setSequence( localIndexes, first, inc, elemsPerPartition );

        dist = DistributionPtr( new GeneralDistribution( globalSize, localIndexes, comm ) );
    }

    ~GeneralDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType elemsPerPartition;
    IndexType globalSize;

    hmemo::HArray<IndexType> localIndexes;

    DistributionPtr dist;
};

BOOST_FIXTURE_TEST_SUITE( GeneralDistributionTest, GeneralDistributionTestConfig )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GeneralDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( generalSizeTest )
{
    BOOST_CHECK( dist->getLocalSize() == elemsPerPartition );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    hmemo::HArray<IndexType> localIndexes;

    IndexType N = 2;

    const IndexType first = rank * N;
    const IndexType inc   = 1;

    utilskernel::HArrayUtils::setSequence( localIndexes, first, inc, N );

    GeneralDistribution genDist1( size * N, localIndexes, comm );
    const GeneralDistribution& genDist2 = genDist1;
    GeneralDistribution genDist3( size * N, localIndexes, comm );

    // general distributions can be compared with each other

    BOOST_CHECK_EQUAL( genDist1, genDist2 );  // pointer equality
    BOOST_CHECK_EQUAL( genDist2, genDist1 );  // pointer equality
    BOOST_CHECK_EQUAL( genDist1, genDist3 );  // same constructor equality
    BOOST_CHECK_EQUAL( genDist3, genDist1 );  // same constructor equality
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getBlockDistributedOwnersTest )
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    const GeneralDistribution* genDist = dynamic_cast<const GeneralDistribution*>( dist.get() );

    BOOST_REQUIRE( genDist != NULL );

    const hmemo::HArray<PartitionId>& localOwners = genDist->getMyBlockDistributedOwners();

    hmemo::HArray<PartitionId> expectedOwners;

    {
        auto wExpected = hmemo::hostWriteOnlyAccess( expectedOwners, elemsPerPartition );

        for ( IndexType i = 0; i < elemsPerPartition; ++i )
        {
            IndexType globalIndex = i + rank * elemsPerPartition;   
            PartitionId owner = globalIndex % size; 
            wExpected[i] = owner;
        }
    }
  
    BOOST_TEST( hostReadAccess( expectedOwners ) == hostReadAccess( localOwners ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistConstructorTest )
{
    IndexType N = 15;

    TestDistributions allDist( N );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        const Distribution& dist = *allDist[i];

        if ( dist.getReduceCommunicator() != *comm )
        {
            // dist is NoDistribution, cannot be redistributed

            continue;
        }

        IndexType nLocal = dist.getLocalSize();

        hmemo::HArray<PartitionId> owners;  // new owners for my local indexes

        {
            hmemo::WriteOnlyAccess<PartitionId> wOwners( owners, nLocal );

            for ( IndexType i = 0; i < nLocal; ++i )
            {
                // choose owner as if it will be a Cyclic(1) distribution

                IndexType globalIndex = dist.local2global( i );
                wOwners[i] = globalIndex % size;
            }
        }

        SCAI_LOG_DEBUG( logger, "redistribute, dist = " << dist << ", owners = " << owners )

        GeneralDistribution gdist( dist, owners );

        CyclicDistribution cyclic( N, 1, comm );

        BOOST_REQUIRE_EQUAL( gdist.getLocalSize(), cyclic.getLocalSize() );

        for ( IndexType i = 0; i < cyclic.getLocalSize(); ++i )
        {
            BOOST_CHECK_EQUAL( cyclic.local2global( i ), gdist.local2global( i ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
