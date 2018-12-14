/**
 * @file GridDistributionTest.cpp
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
 * @brief Specific tests for the class GridDistribution.
 * @author Thomas Brandes
 * @date 30.01.2017
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/GridDistribution.hpp>
#include <scai/utilskernel.hpp>

using namespace scai;
using namespace dmemo;
using common::Grid;
using common::Grid1D;
using common::Grid2D;
using common::Grid3D;

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
        dist = DistributionPtr( new GridDistribution( Grid1D( blockSize * size ), comm, Grid1D( size ) ) );
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
    bdist = Distribution::getDistributionPtr( "GRID", comm, globalSize );
    BOOST_CHECK_EQUAL( bdist->getGlobalSize(), globalSize );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    // do not use the global comm variable here

    CommunicatorPtr comm;

    const IndexType globalSize = 17;

    Grid2D globalGrid( globalSize, globalSize );

    // constructor must fail for empty communicator

    BOOST_CHECK_THROW(
    {
        GridDistribution bdist( globalGrid, comm, Grid1D( 1 ) );
    },
    common::Exception );

    comm = Communicator::getCommunicatorPtr();

    GridDistribution bdist( globalGrid, comm );

    BOOST_CHECK_EQUAL( bdist.getGlobalSize(), globalSize * globalSize );

    BOOST_CHECK_EQUAL( globalGrid, bdist.getGlobalGrid() );

    const Grid& localGrid = bdist.getLocalGrid();

    // local grid must be subgrid

    BOOST_CHECK_EQUAL( globalGrid.nDims(), localGrid.nDims() );

    for ( IndexType idim = 0; idim < globalGrid.nDims(); ++ idim )
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

    Grid2D globalGrid( N1, N2 );
    Grid2D procGrid( 1, 1 );

    GridDistribution myGridDistribution( globalGrid, comm, procGrid );

    const Grid& localGrid = myGridDistribution.getLocalGrid();

    if ( comm->getRank() == 0 )
    {
        BOOST_CHECK_EQUAL( localGrid, globalGrid );
    }
    else
    {
        const IndexType zero = 0;
        BOOST_CHECK_EQUAL( localGrid, Grid2D( zero, zero ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructor2Test )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType globalSize = 17;

    Grid3D globalGrid( globalSize, globalSize, globalSize );
    Grid3D procGrid( comm->getSize(), 1, 1 );

    GridDistribution myGridDistribution( globalGrid, comm, procGrid );

    const Grid& localGrid = myGridDistribution.getLocalGrid();

    BOOST_CHECK_EQUAL( globalGrid.nDims(), localGrid.nDims() );

    if ( comm->getSize() == 1 )
    {
        BOOST_CHECK_EQUAL( localGrid.size( 0 ), globalGrid.size( 0 ) );
    }
    else
    {
        // first dimension must be distributed

        BOOST_CHECK( localGrid.size( 0 ) < globalGrid.size( 0 ) );
    }

    for ( IndexType idim = 1; idim < globalGrid.nDims(); ++ idim )
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
        GridDistribution small( Grid1D( n ), comm, Grid1D( comm->getSize() ) );

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
    using namespace hmemo;

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
    using namespace hmemo;

    HArray<IndexType> myIndexes1;
    HArray<IndexType> myIndexes2;

    dist->getOwnedIndexes( myIndexes1 );                 // call it for block distribution
    dist->Distribution::getOwnedIndexes( myIndexes2 );   // call if from base class

    IndexType nLocal = dist->getLocalSize();

    BOOST_REQUIRE_EQUAL( nLocal, myIndexes1.size() );
    BOOST_REQUIRE_EQUAL( nLocal, myIndexes2.size() );

    BOOST_TEST( hostReadAccess( myIndexes1 ) == hostReadAccess( myIndexes2 ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    Grid2D globalGrid( 10, 10 );
    Grid2D procGrid( 1, comm->getSize() );

    DistributionPtr griddist1( new GridDistribution( globalGrid, comm, procGrid ) );
    DistributionPtr griddist2( griddist1 );
    DistributionPtr griddist3( new GridDistribution( globalGrid, comm, procGrid ) );
    BOOST_CHECK( ( *griddist1 ).isEqual( *griddist2 ) );
    BOOST_CHECK( ( *griddist1 ).isEqual( *griddist3 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
