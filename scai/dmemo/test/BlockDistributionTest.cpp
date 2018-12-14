/**
 * @file BlockDistributionTest.cpp
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
 * @brief Contains the implementation of the class BlockDistributionTest.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/utilskernel.hpp>

using namespace scai;
using namespace dmemo;

using hmemo::HArray;

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
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    const IndexType globalSize = 5;

    DistributionPtr bdist = BlockDistribution::create( DistributionArguments( comm, globalSize, NULL, 1.0 ) );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), globalSize );
    bdist = Distribution::getDistributionPtr( "BLOCK", comm, globalSize );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), globalSize );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    CommunicatorPtr comm;
    const IndexType globalSize = 17;

    BOOST_CHECK_THROW(
    {
        BlockDistribution bdist( globalSize, comm );
    },
    common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( blockSizeTest )
{
    // Test smaller sizes
    for ( IndexType n = 1; n <= static_cast<IndexType>( size ); n++ )
    {
        BlockDistribution small( n, comm );

        // only the first n partitions have one element

        if ( static_cast<IndexType>( rank ) < n )
        {
            BOOST_CHECK_EQUAL( small.getLocalSize(), IndexType( 1 ) );
        }
        else
        {
            BOOST_CHECK_EQUAL( small.getLocalSize(), IndexType( 0 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( blockComputeOwnersTest )
{
    using namespace utilskernel;

    HArray<IndexType> indexes;
    HArrayUtils::setOrder( indexes, size * blockSize );

    HArray<PartitionId> owners;
    dist->computeOwners( owners, indexes );

    HArray<PartitionId> owners1;

    BOOST_CHECK_EQUAL( owners.size(), indexes.size() );

    hmemo::ReadAccess<PartitionId> rOwners( owners );

    for ( PartitionId p = 0; p < size; ++p )
    {
        for ( IndexType i = 0; i < blockSize; ++i )
        {
            IndexType pos = p * blockSize + i;

            BOOST_CHECK_EQUAL( p, owners[pos] );

            if ( p == rank )
            {
                BOOST_CHECK( dist->isLocal( pos ) );
                BOOST_CHECK_EQUAL( dist->global2Local( pos ), i );
            }
            else
            {
                BOOST_CHECK( !dist->isLocal( pos ) );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ownedIndexesTest )
{
    SCAI_LOG_INFO( logger, "ownedIndexesTest for dist = " << *dist )

    using namespace utilskernel;

    HArray<IndexType> myIndexes1;
    HArray<IndexType> myIndexes2;

    dist->getOwnedIndexes( myIndexes1 );                 // call it for block distribution
    dist->Distribution::getOwnedIndexes( myIndexes2 );   // call if from base class

    IndexType nLocal = dist->getLocalSize();

    BOOST_REQUIRE_EQUAL( nLocal, myIndexes1.size() );
    BOOST_REQUIRE_EQUAL( nLocal, myIndexes2.size() );

    hmemo::ReadAccess<IndexType> rIndexes1( myIndexes1 );
    hmemo::ReadAccess<IndexType> rIndexes2( myIndexes2 );

    for ( IndexType i = 0; i < nLocal; ++i )
    {
        BOOST_CHECK_EQUAL( rIndexes1[i], rIndexes2[i] );
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
