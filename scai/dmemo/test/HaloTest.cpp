/**
 * @file HaloTest.cpp
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
 * @brief Test routines for the class Halo.
 * @author Thomas Brandes
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/logging.hpp>

#include <scai/dmemo/Halo.hpp>
#include <scai/dmemo/HaloBuilder.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/BlockDistribution.hpp>

#include <scai/common/test/TestMacros.hpp>

using namespace scai;
using namespace dmemo;
using namespace hmemo;

/* --------------------------------------------------------------------- */

/** Fixture that is called with each test of this test suite. */

struct HaloTestConfig
{
    HaloTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();

        rank = comm->getRank();
        size = comm->getSize();

        dist.reset( new BlockDistribution( size, comm ) );

        leftNeighbor = comm->getNeighbor( -1 );
        rightNeighbor = comm->getNeighbor( 1 );

        buildRequiredIndexes();

        {
            HArrayRef<IndexType> arrRequiredIndexes( requiredIndexes );
            HaloBuilder::build( *dist, arrRequiredIndexes, halo );
        }
    }

    ~HaloTestConfig()
    {
    }

    CommunicatorPtr comm;
    DistributionPtr dist;
    std::vector<IndexType> requiredIndexes;
    PartitionId leftNeighbor;
    PartitionId rightNeighbor;
    PartitionId rank;
    PartitionId size;
    Halo halo;

    void buildRequiredIndexes()
    {
        // Each processor requires values from left and right neighbor

        if ( !dist->isLocal( leftNeighbor ) )
        {
            requiredIndexes.push_back( leftNeighbor );
        }

        if ( rightNeighbor != leftNeighbor && !dist->isLocal( rightNeighbor ) )
        {
            requiredIndexes.push_back( rightNeighbor );
        }
    }
};

/* --------------------------------------------------------------------- */

BOOST_FIXTURE_TEST_SUITE( HaloTest, HaloTestConfig )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HaloTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildHaloTest )
{
    const IndexType noReqIndexes = static_cast<IndexType>( requiredIndexes.size() );

    const Halo& haloRef = halo;
    const CommunicationPlan& requiredPlan = haloRef.getRequiredPlan();
    const CommunicationPlan& providesPlan = haloRef.getProvidesPlan();

    // check for a correct provide plan

    IndexType offsetCheck = 0;

    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        IndexType n = requiredPlan[p].quantity;
        BOOST_CHECK_EQUAL( ( IndexType ) 1, n );
        PartitionId neighbor = requiredPlan[p].partitionId;
        BOOST_CHECK( neighbor == leftNeighbor || neighbor == rightNeighbor );
        BOOST_CHECK_EQUAL( requiredPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( noReqIndexes, requiredPlan.totalQuantity() );

    // check for a correct provide plan

    offsetCheck = 0;
    PartitionId nProvides = providesPlan.size();

    for ( PartitionId p = 0; p < nProvides; ++p )
    {
        IndexType n = providesPlan[p].quantity;
        BOOST_CHECK_EQUAL( n, static_cast<IndexType>( 1 ) );
        PartitionId neighbor = providesPlan[p].partitionId;
        BOOST_CHECK( neighbor == leftNeighbor || neighbor == rightNeighbor );
        BOOST_CHECK_EQUAL( providesPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    // we have a symmetric halo exchange
    BOOST_CHECK_EQUAL( noReqIndexes, providesPlan.totalQuantity() );

    const ReadAccess<IndexType> providesIndexes( haloRef.getProvidesIndexes() );

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        const IndexType* indexes = providesIndexes + providesPlan[p].offset;
        IndexType expectedLocalIndex = rank;
        BOOST_CHECK_EQUAL( expectedLocalIndex, dist->local2global( indexes[0] ) );
    }

    BOOST_CHECK_EQUAL( noReqIndexes, halo.getHaloSize() );
    IndexType nIndexes = static_cast<IndexType>( requiredIndexes.size() );

    for ( IndexType i = 0; i < nIndexes; ++i )
    {
        const IndexType haloIndex = halo.global2halo( requiredIndexes[i] );
        BOOST_CHECK( common::Utils::validIndex( haloIndex, halo.getHaloSize() ) );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyHaloTest )
{
    const Halo& halo1 = halo;
    Halo halo2( halo1 );  // COPY constructor

    const CommunicationPlan& providesPlan1 = halo1.getProvidesPlan();
    const CommunicationPlan& providesPlan2 = halo2.getProvidesPlan();

    BOOST_REQUIRE_EQUAL( providesPlan1.size(), providesPlan2.size() );
    BOOST_REQUIRE_EQUAL( providesPlan1.totalQuantity(), providesPlan2.totalQuantity() );

    // check equality of both provides plans

    PartitionId nProvides = providesPlan1.size();

    for ( PartitionId p = 0; p < nProvides; ++p )
    {
        BOOST_CHECK_EQUAL( providesPlan1[p].quantity, providesPlan2[p].quantity );
        BOOST_CHECK_EQUAL( providesPlan1[p].partitionId, providesPlan2[p].partitionId );
    };

    IndexType nIndexes = providesPlan1.totalQuantity();

    BOOST_CHECK_EQUAL( halo1.getProvidesIndexes().size(), halo2.getProvidesIndexes().size() );

    {
        const ReadAccess<IndexType> providesIndexes1( halo1.getProvidesIndexes() );
        const ReadAccess<IndexType> providesIndexes2( halo2.getProvidesIndexes() );

        for ( IndexType i = 0; i < nIndexes; ++i )
        {
            BOOST_CHECK_EQUAL( providesIndexes1[i], providesIndexes2[i] );
        }

        BOOST_CHECK_THROW(
        {
            halo2.clear();
        }, common::Exception );
    }

    // ReadAccesses are released, so clear is safe now

    halo.clear();  // -> halo1.clear()

    BOOST_CHECK_EQUAL( IndexType( 0 ), halo1.getProvidesIndexes().size() );
    BOOST_CHECK_EQUAL( nIndexes, halo2.getProvidesIndexes().size() );

    halo2.purge();

    BOOST_CHECK_EQUAL( IndexType( 0 ), halo2.getProvidesIndexes().size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( updateHaloTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();
    SCAI_LOG_INFO( logger, "updateHaloTest<" << common::getScalarType<ValueType>() << ">" );
    const IndexType factor = 4;
    const IndexType vectorSize = factor * size;
    BlockDistribution distribution( vectorSize, comm );
    std::vector<IndexType> requiredIndexes;

    for ( IndexType i = 0; i < factor; ++i )
    {
        const IndexType requiredIndex = ( ( rank + 1 ) * factor + i ) % vectorSize;

        if ( distribution.isLocal( requiredIndex ) )
        {
            continue;
        }

        requiredIndexes.push_back( requiredIndex );
    }

    SCAI_LOG_INFO( logger, "build the Halo" );

    Halo halo;
    {
        HArrayRef<IndexType> arrRequiredIndexes( requiredIndexes );
        HaloBuilder::build( distribution, arrRequiredIndexes, halo );
    }

    SCAI_LOG_INFO( logger, "halo is now available: " << halo );
    HArray<ValueType> localData;
    {
        WriteOnlyAccess<ValueType> localDataAccess( localData, distribution.getLocalSize() );

        for ( IndexType i = 0; i < localData.size(); ++i )
        {
            localDataAccess[i] = static_cast<ValueType>( distribution.local2global( i ) );
        }
    }
    SCAI_LOG_INFO( logger, "update halo data by communicator" );
    HArray<ValueType> haloData;
    comm->updateHalo( haloData, localData, halo );
    BOOST_CHECK_EQUAL( static_cast<IndexType>( requiredIndexes.size() ), haloData.size() );
    {
        ReadAccess<ValueType> haloDataAccess( haloData );

        for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexes.size() ); ++i )
        {
            ValueType expectedValue = static_cast<ValueType>( requiredIndexes[i] );
            BOOST_CHECK_EQUAL( expectedValue, haloDataAccess[i] );
        }
    }
    requiredIndexes.clear();

    for ( IndexType i = 0; i < vectorSize; ++i )
    {
        if ( distribution.isLocal( i ) || ( i + rank ) % 2 == 0 )
        {
            continue;
        }

        requiredIndexes.push_back( i );
    }

    {
        HArrayRef<IndexType> arrRequiredIndexes( requiredIndexes );
        HaloBuilder::build( distribution, arrRequiredIndexes, halo );
    }

    {
        common::unique_ptr<tasking::SyncToken> token( comm->updateHaloAsync( haloData, localData, halo ) );
    }

    BOOST_CHECK_EQUAL( static_cast<IndexType>( requiredIndexes.size() ), haloData.size() );
    {
        ReadAccess<ValueType> haloDataAccess( haloData );

        for ( IndexType i = 0; i < static_cast<IndexType>( requiredIndexes.size() ); ++i )
        {
            ValueType expectedValue = static_cast<ValueType>( requiredIndexes[i] );
            BOOST_CHECK_EQUAL( expectedValue, haloDataAccess[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
