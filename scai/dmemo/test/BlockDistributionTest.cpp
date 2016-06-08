/**
 * @file BlockDistributionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief Contains the implementation of the class BlockDistributionTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/BlockDistribution.hpp>

using namespace scai::dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct BlockDistributionTestConfig
{
    BlockDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
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

/* --------------------------------------------------------------------- */

BOOST_FIXTURE_TEST_SUITE( BlockDistributionTest, BlockDistributionTestConfig );

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.BlockDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    // verify that a BlockDistribution can be created via the factory

    IndexType globalSize = 17;

    DistributionPtr bdist( Distribution::getDistributionPtr( "BLOCK", comm, globalSize, 1.0 ) );

    BOOST_REQUIRE( bdist );

    SCAI_LOG_INFO( logger, "created by factory: " << *bdist )

    // verify by dynamic cast that it is really a BlockDistribution

    const BlockDistribution* b = dynamic_cast<const BlockDistribution*>( bdist.get() );

    BOOST_CHECK( b );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    DistributionPtr bdist ( BlockDistribution::create( DistributionArguments( comm, 1, NULL, 1.0 ) ) );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), 1 );
    bdist.reset( Distribution::getDistributionPtr( "BLOCK", comm, 1 ) );
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
