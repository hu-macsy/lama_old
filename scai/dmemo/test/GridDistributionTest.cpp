/**
 * @file GridDistributionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2017
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Specific tests for the class GridDistribution.
 * @author Thomas Brandes
 * @date 30.01.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/utilskernel/LArray.hpp>

using namespace scai;
using namespace dmemo;
using common::Grid;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct GridDistributionTestConfig
{
    GridDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
        rank = comm->getRank();
        size = comm->getSize();
        blockSize = 17;
        dist = DistributionPtr( new GridDistribution( Grid( blockSize * size ), comm, Grid( size ) ) );
    }

    ~GridDistributionTestConfig()
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

BOOST_FIXTURE_TEST_SUITE( GridDistributionTest, GridDistributionTestConfig );

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GridDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    // verify that a GridDistribution can be created via the factory
    IndexType globalSize = 17;
    DistributionPtr bdist( Distribution::getDistributionPtr( "GRID", comm, globalSize, 1.0 ) );
    BOOST_REQUIRE( bdist );
    SCAI_LOG_INFO( logger, "created by factory: " << *bdist )
    // verify by dynamic cast that it is really a GridDistribution
    const GridDistribution* b = dynamic_cast<const GridDistribution*>( bdist.get() );
    BOOST_CHECK( b );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    const IndexType globalSize = 5;

    DistributionPtr bdist ( GridDistribution::create( DistributionArguments( comm, globalSize, NULL, 1.0 ) ) );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), globalSize );
    bdist.reset( Distribution::getDistributionPtr( "GRID", comm, globalSize ) );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), globalSize );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    // do not use the global comm variable here

    CommunicatorPtr comm;

    const IndexType globalSize = 17;

    Grid globalGrid( globalSize, globalSize );

    // constructor must fail for empty communicator

    BOOST_CHECK_THROW(
    {
        GridDistribution bdist( globalGrid, comm, Grid( 1 ) );
    },
    common::Exception );

    comm = Communicator::getCommunicatorPtr();

    GridDistribution bdist( globalGrid, comm );

    BOOST_CHECK_EQUAL( bdist.getGlobalSize(), globalSize * globalSize );

    BOOST_CHECK_EQUAL( globalGrid, bdist.getGlobalGrid() );

    const Grid& localGrid = bdist.getLocalGrid();

    // local grid must be subgrid

    BOOST_CHECK_EQUAL( globalGrid.ndims(), localGrid.ndims() );

    for ( IndexType idim = 0; idim < globalGrid.ndims(); ++ idim )
    {
        BOOST_CHECK( localGrid.size( idim ) <= globalGrid.size( idim ) );
    }

    if ( comm->getSize() == 1 )
    {
        // only one processor available, local and global grid must be same

        BOOST_CHECK_EQUAL( globalGrid, localGrid );
    }
    else
    {
        // grid is distributed in any case

        BOOST_CHECK( globalGrid != localGrid );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructor1Test )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType N1 = 1;
    const IndexType N2 = 2;

    Grid globalGrid( N1, N2 );
    Grid procGrid( 1, 1 );

    GridDistribution myGridDistribution( globalGrid, comm, procGrid );

    const Grid& localGrid = myGridDistribution.getLocalGrid();

    if ( comm->getRank() == 0 )
    {
        BOOST_CHECK_EQUAL( localGrid, globalGrid );
    }
    else
    {
        BOOST_CHECK_EQUAL( localGrid, Grid( 0, 0 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructor2Test )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType globalSize = 17;

    Grid globalGrid( globalSize, globalSize, globalSize );
    Grid procGrid( comm->getSize(), 1, 1 );

    GridDistribution myGridDistribution( globalGrid, comm, procGrid );

    const Grid& localGrid = myGridDistribution.getLocalGrid();

    BOOST_CHECK_EQUAL( globalGrid.ndims(), localGrid.ndims() );

    if ( comm->getSize() == 1 )
    {
        BOOST_CHECK_EQUAL( localGrid.size( 0 ), globalGrid.size( 0 ) );
    }
    else
    {
        // first dimension must be distributed

        BOOST_CHECK( localGrid.size( 0 ) < globalGrid.size( 0 ) );
    }

    for ( IndexType idim = 1; idim < globalGrid.ndims(); ++ idim )
    {
        // all other dimensions are not distributed at all

        BOOST_CHECK_EQUAL( localGrid.size( idim ), globalGrid.size( idim ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( blockSizeTest )
{
    // Test smaller sizes
    for ( IndexType n = 1; n <= static_cast<IndexType>( size ); n++ )
    {
        GridDistribution small( Grid( n ), comm, Grid( comm->getSize() ) );

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

    LArray<IndexType> indexes;
    HArrayUtils::setOrder( indexes, size * blockSize );

    LArray<PartitionId> owners;
    dist->computeOwners( owners, indexes );

    LArray<PartitionId> owners1;

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
                BOOST_CHECK_EQUAL( dist->global2local( pos ), i );
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

    LArray<IndexType> myIndexes1;
    LArray<IndexType> myIndexes2;

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
    Grid globalGrid( 10, 10 );
    Grid procGrid( 1, comm->getSize() );

    DistributionPtr griddist1( new GridDistribution( globalGrid, comm, procGrid ) );
    DistributionPtr griddist2( griddist1 );
    DistributionPtr griddist3( new GridDistribution( globalGrid, comm, procGrid ) );
    BOOST_CHECK( ( *griddist1 ).isEqual( *griddist2 ) );
    BOOST_CHECK( ( *griddist1 ).isEqual( *griddist3 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();