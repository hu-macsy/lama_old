/**
 * @file HaloPlanTest.cpp
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
 * @brief Test routines for the class HaloPlan
 * @author Thomas Brandes, Andreas Longva
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

#include <scai/dmemo/HaloPlan.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/dmemo/test/TestMacros.hpp>

#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>

#include <map>

using namespace scai;
using namespace dmemo;
using namespace hmemo;

using namespace boost::test_tools;

/* --------------------------------------------------------------------- */

/** Fixture that is called with each test of this test suite. */

struct HaloPlanTestConfig
{
    HaloPlanTestConfig()
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
            haloPlan = haloPlanByRequiredIndexes( arrRequiredIndexes, *dist );
        }
    }

    ~HaloPlanTestConfig()
    {
    }

    CommunicatorPtr comm;
    DistributionPtr dist;
    std::vector<IndexType> requiredIndexes;
    PartitionId leftNeighbor;
    PartitionId rightNeighbor;
    PartitionId rank;
    PartitionId size;
    HaloPlan haloPlan;

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

BOOST_FIXTURE_TEST_SUITE( HaloPlanTest, HaloPlanTestConfig )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HaloPlanTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildHaloPlanTest )
{
    const IndexType noReqIndexes = static_cast<IndexType>( requiredIndexes.size() );

    const HaloPlan& haloRef = haloPlan;
    const CommunicationPlan& requiredPlan = haloRef.getHaloCommunicationPlan();
    const CommunicationPlan& providesPlan = haloRef.getLocalCommunicationPlan();

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

    /**
    const ReadAccess<IndexType> providesIndexes( haloRef.getProvidesIndexes() );

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        providesPlan.getInfo( quantity, offset
        const IndexType* indexes = providesIndexes + providesPlan[p].offset;
        IndexType expectedLocalIndex = rank;
        BOOST_CHECK_EQUAL( expectedLocalIndex, dist->local2global( indexes[0] ) );
    }
    */

    BOOST_CHECK_EQUAL( noReqIndexes, haloPlan.getHaloSize() );
    IndexType numIndexes = static_cast<IndexType>( requiredIndexes.size() );

    for ( IndexType i = 0; i < numIndexes; ++i )
    {
        const IndexType haloIndex = haloPlan.global2halo( requiredIndexes[i] );
        SCAI_LOG_INFO( logger, "check valid halo index = " << haloIndex << " was global " << requiredIndexes[i] 
                                << " for haloSize = " << haloPlan.getHaloSize() )
        BOOST_CHECK( common::Utils::validIndex( haloIndex, haloPlan.getHaloSize() ) );
    }

    for ( IndexType i = 0; i < 10; ++i )
    {
        const IndexType haloIndex = haloPlan.global2halo( i );

        if ( std::find( requiredIndexes.begin(), requiredIndexes.end(), i ) == requiredIndexes.end() )
        {
            BOOST_CHECK_EQUAL( haloIndex, invalidIndex );
        }
    }
}

// Helper struct to simplify building expected values
// on multiple processors
struct HaloExpectedResult
{
    std::vector<IndexType> providedIndexes;
    std::vector<IndexType> requiredIndexes;
    std::vector<IndexType> providedQuantities;
    std::vector<IndexType> requiredQuantities;
    std::map<IndexType, IndexType> global2halo;
};

void checkHaloAgainstExpected( const HaloPlan& halo, const HaloExpectedResult& expected )
{
    const auto expectedProvidedPlan = CommunicationPlan( expected.providedQuantities );
    const auto expectedRequiredPlan = CommunicationPlan( expected.requiredQuantities );

    CHECK_COMMUNICATION_PLANS_EQUAL( expectedProvidedPlan, halo.getLocalCommunicationPlan() );
    CHECK_COMMUNICATION_PLANS_EQUAL( expectedRequiredPlan, halo.getHaloCommunicationPlan() );

    BOOST_TEST( hostReadAccess( halo.getProvidesIndexes() ) == expected.providedIndexes, boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( halo.getHalo2GlobalIndexes() ) == expected.requiredIndexes, boost::test_tools::per_element() );

    BOOST_TEST_CONTEXT( " for global2halo comparison" )
    {
        for ( auto keyValue : expected.global2halo )
        {
            const auto globalIndex = keyValue.first;
            const auto haloIndex = keyValue.second;
            BOOST_TEST ( halo.global2halo( globalIndex ) == haloIndex );
        }
    }
};

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyHaloPlanTest )
{
    const HaloPlan& halo1 = haloPlan;
    HaloPlan halo2( halo1 );  // COPY constructor

    const CommunicationPlan& providesPlan1 = halo1.getLocalCommunicationPlan();
    const CommunicationPlan& providesPlan2 = halo2.getLocalCommunicationPlan();

    BOOST_REQUIRE_EQUAL( providesPlan1.size(), providesPlan2.size() );
    BOOST_REQUIRE_EQUAL( providesPlan1.totalQuantity(), providesPlan2.totalQuantity() );

    // check equality of both provides plans

    PartitionId nProvides = providesPlan1.size();

    for ( PartitionId p = 0; p < nProvides; ++p )
    {
        BOOST_CHECK_EQUAL( providesPlan1[p].quantity, providesPlan2[p].quantity );
        BOOST_CHECK_EQUAL( providesPlan1[p].partitionId, providesPlan2[p].partitionId );
    };

    IndexType numIndexes = providesPlan1.totalQuantity();

    BOOST_CHECK_EQUAL( halo1.getProvidesIndexes().size(), halo2.getProvidesIndexes().size() );

    {
        const ReadAccess<IndexType> providesIndexes1( halo1.getProvidesIndexes() );
        const ReadAccess<IndexType> providesIndexes2( halo2.getProvidesIndexes() );

        for ( IndexType i = 0; i < numIndexes; ++i )
        {
            BOOST_CHECK_EQUAL( providesIndexes1[i], providesIndexes2[i] );
        }

        BOOST_CHECK_THROW(
        {
            halo2.clear();
        }, common::Exception );
    }

    // ReadAccesses are released, so clear is safe now

    haloPlan.clear();  // -> halo1.clear()

    BOOST_CHECK_EQUAL( IndexType( 0 ), halo1.getProvidesIndexes().size() );
    BOOST_CHECK_EQUAL( numIndexes, halo2.getProvidesIndexes().size() );

    halo2.purge();

    BOOST_CHECK_EQUAL( IndexType( 0 ), halo2.getProvidesIndexes().size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE_TEMPLATE( updateHaloPlanTest, ValueType, scai_numeric_test_types )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
    BOOST_REQUIRE( comm );
    IndexType rank = comm->getRank();
    IndexType size = comm->getSize();
    SCAI_LOG_INFO( logger, "updateHaloPlanTest<" << common::getScalarType<ValueType>() << ">" );
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

    HaloPlan halo;
    {
        HArrayRef<IndexType> arrRequiredIndexes( requiredIndexes );
        halo = haloPlanByRequiredIndexes( arrRequiredIndexes, distribution );
    }

    SCAI_LOG_INFO( logger, "halo is now available: " << halo );
    HArray<ValueType> localData;
    {
        WriteOnlyAccess<ValueType> localDataAccess( localData, distribution.getLocalSize() );

        for ( IndexType i = 0; i < localData.size(); ++i )
        {
            localDataAccess[i] = static_cast<ValueType>( distribution.local2Global( i ) );
        }
    }
    SCAI_LOG_INFO( logger, "update halo data by communicator" );
    HArray<ValueType> haloData;
    halo.updateHalo( haloData, localData, *comm );
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
        haloPlan = haloPlanByRequiredIndexes( arrRequiredIndexes, distribution );
    }

    {
        std::unique_ptr<tasking::SyncToken> token( haloPlan.updateHaloAsync( haloData, localData, *comm ) );
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

BOOST_AUTO_TEST_CASE( exampleTest )
{

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
 
    if ( comm->getSize() != 3 ) 
    {
        return;
    }

    IndexType localSize = 0;

    HArray<IndexType> requiredIndexes;

    HArray<IndexType> expHalo2GlobalIndexes;
    HArray<IndexType> expProvidesIndexes;

    HArray<IndexType> expProvidesSizes;
    HArray<IndexType> expHaloSizes;

    // Communication-Matrix
    //
    //           0     1    2
    //      0    -     2    2        -   3,4  6,7    -  0,1  0,1
    //      1    1     -    2        0    -   7,9    0    -  1,3
    //      2    1     2    -        2   3,5   -     2  0,2   - 
    if ( comm->getRank() == 0 )
    {
        localSize = 3;
        requiredIndexes = HArray<IndexType>( { 3, 6, 4, 7 } );

        expHaloSizes = HArray<IndexType>( { 0, 2, 2  } );
        expHalo2GlobalIndexes = HArray<IndexType>( { 3, 4, 6, 7 } );

        expProvidesSizes = HArray<IndexType>( { 0, 1, 1  } );
        expProvidesIndexes = HArray<IndexType>( { 0, 2 } );
    }
    else if ( comm->getRank() == 1 )
    {
        localSize = 3;
        requiredIndexes = HArray<IndexType>( { 0, 7, 9 } );

        expHaloSizes = HArray<IndexType>( { 1, 0, 2  } );
        expHalo2GlobalIndexes = HArray<IndexType>( { 0, 7, 9 } );
        expProvidesSizes = HArray<IndexType>( { 2, 0, 2 } );
        expProvidesIndexes = HArray<IndexType>( { 0, 1, 0, 2 } );
    }
    else if ( comm->getRank() == 2 )
    {
        localSize = 4;
        requiredIndexes = HArray<IndexType>( { 2, 3, 5 } );

        expHaloSizes = HArray<IndexType>( { 1, 2, 0  } );
        expHalo2GlobalIndexes = HArray<IndexType>( { 2, 3, 5 } );

        expProvidesSizes = HArray<IndexType>( { 2, 2, 0 } );
        expProvidesIndexes = HArray<IndexType>( { 0, 1, 1, 3 } );
    }

    auto blockDist = genBlockDistribution( localSize, comm );

    auto haloPlan = haloPlanByRequiredIndexes( requiredIndexes, *blockDist );
    
    BOOST_TEST( hostReadAccess( expHalo2GlobalIndexes ) == hostReadAccess( haloPlan.getHalo2GlobalIndexes()), per_element() );
    BOOST_TEST( hostReadAccess( expProvidesIndexes ) == hostReadAccess( haloPlan.getProvidesIndexes()), per_element() );

    CommunicationPlan expRequiredPlan( hostReadAccess( expHaloSizes ) );
    CommunicationPlan expProvidesPlan( hostReadAccess( expProvidesSizes ) );

    CHECK_COMMUNICATION_PLANS_EQUAL( expRequiredPlan, haloPlan.getHaloCommunicationPlan() );
    CHECK_COMMUNICATION_PLANS_EQUAL( expProvidesPlan, haloPlan.getLocalCommunicationPlan() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
