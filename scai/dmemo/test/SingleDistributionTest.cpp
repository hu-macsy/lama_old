/**
 * @file SingleDistributionTest.cpp
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
 * @brief Contains the implementation of the class SingleDistributionTest.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/SingleDistribution.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct SingleDistributionTestConfig
{
    SingleDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
        rank = comm->getRank();
        size = comm->getSize();
        blockSize = 17;
        dist = DistributionPtr( new SingleDistribution( blockSize * size, comm, 0 ) );
    }

    ~SingleDistributionTestConfig()
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

BOOST_FIXTURE_TEST_SUITE( SingleDistributionTest, SingleDistributionTestConfig );

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.SingleDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( factoryTest )
{
    // verify that a SingleDistribution can be created via the factory
    IndexType globalSize = 17;
    DistributionPtr bdist( Distribution::getDistributionPtr( "SINGLE", comm, globalSize, 1.0 ) );
    BOOST_REQUIRE( bdist );
    SCAI_LOG_INFO( logger, "created by factory: " << *bdist )
    // verify by dynamic cast that it is really a SingleDistribution
    const SingleDistribution* b = dynamic_cast<const SingleDistribution*>( bdist.get() );
    BOOST_CHECK( b );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( createTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    const IndexType globalSize = 5;

    DistributionPtr sdist ( SingleDistribution::create( DistributionArguments( comm, globalSize, NULL, 1.0 ) ) );
    BOOST_CHECK_EQUAL( sdist->getGlobalSize(), globalSize );
    sdist.reset( Distribution::getDistributionPtr( "SINGLE", comm, globalSize ) );
    BOOST_CHECK_EQUAL( sdist->getGlobalSize(), globalSize );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    CommunicatorPtr comm;
    const IndexType globalSize = 17;

    BOOST_CHECK_THROW(
    {
        SingleDistribution bdist( globalSize, comm, 0 );
    },
    common::Exception );

    comm = Communicator::getCommunicatorPtr();
    BOOST_CHECK_THROW(
    {
        SingleDistribution bdist( globalSize, comm, comm->getSize() );
    },
    common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( localSizeTest )
{
    // Test smaller sizes
    for ( IndexType n = 1; n <= static_cast<IndexType>( size ); n++ )
    {
        PartitionId owner = 0;

        SingleDistribution small( n, comm, owner );

        // only the first n partitions have one element

        if ( comm->getRank() == owner )
        {
            BOOST_CHECK_EQUAL( small.getLocalSize(), n );
        }
        else
        {
            BOOST_CHECK_EQUAL( small.getLocalSize(), IndexType( 0 ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ownedIndexesTest )
{
    SCAI_LOG_INFO( logger, "ownedIndexesTest for dist = " << *dist )

    using namespace hmemo;

    HArray<IndexType> myIndexes1;
    HArray<IndexType> myIndexes2;

    dist->getOwnedIndexes( myIndexes1 );                 // call it for single distribution
    dist->Distribution::getOwnedIndexes( myIndexes2 );   // call if from base class

    IndexType nLocal = dist->getLocalSize();

    BOOST_REQUIRE_EQUAL( nLocal, myIndexes1.size() );
    BOOST_REQUIRE_EQUAL( nLocal, myIndexes2.size() );

    BOOST_TEST( hostReadAccess( myIndexes1 ) == hostReadAccess( myIndexes2 ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    const IndexType n = 5;

    DistributionPtr singledist1( new SingleDistribution( n, comm, 0 ) );
    DistributionPtr singledist2( singledist1 );
    DistributionPtr singledist3( new SingleDistribution( n, comm, 0 ) );
    DistributionPtr singledist4( new SingleDistribution( n - 1, comm, 0 ) );
    BOOST_CHECK( ( *singledist1 ).isEqual( *singledist2 ) );
    BOOST_CHECK( ( *singledist1 ).isEqual( *singledist3 ) );
    BOOST_CHECK( !( *singledist1 ).isEqual( *singledist4 ) );

    if ( comm->getSize() > 1 )
    {
        DistributionPtr singledist5( new SingleDistribution( n, comm, 1 ) );
        BOOST_CHECK( !( *singledist1 ).isEqual( *singledist5 ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
