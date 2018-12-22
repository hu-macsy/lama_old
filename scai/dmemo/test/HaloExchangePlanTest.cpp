/**
 * @file HaloExchangePlanTest.cpp
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
 * @brief Test routines for the class HaloExchangePlan
 * @author Thomas Brandes, Andreas Longva
 * @date 24.07.2016
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>

#include <scai/common/SCAITypes.hpp>

#include <scai/logging.hpp>

#include <scai/dmemo/HaloExchangePlan.hpp>
#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/dmemo/test/TestMacros.hpp>

#include <scai/hmemo/HostReadAccess.hpp>
#include <scai/hmemo/HostWriteAccess.hpp>

#include <scai/utilskernel/HArrayUtils.hpp>

#include <map>

using namespace scai;
using namespace dmemo;
using namespace hmemo;

using namespace boost::test_tools;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( HaloExchangePlanTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.HaloExchangePlanTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( exampleTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();
 
    // This test runs only on exact 3 processors

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

    auto haloExchangePlan = haloExchangePlanByRequiredIndexes( requiredIndexes, *blockDist );
    
    BOOST_TEST( hostReadAccess( expHalo2GlobalIndexes ) == hostReadAccess( haloExchangePlan.getHalo2GlobalIndexes()), per_element() );
    BOOST_TEST( hostReadAccess( expProvidesIndexes ) == hostReadAccess( haloExchangePlan.getLocalIndexes()), per_element() );

    CommunicationPlan expRequiredPlan( hostReadAccess( expHaloSizes ) );
    CommunicationPlan expProvidesPlan( hostReadAccess( expProvidesSizes ) );

    CHECK_COMMUNICATION_PLANS_EQUAL( expRequiredPlan, haloExchangePlan.getHaloCommunicationPlan() );
    CHECK_COMMUNICATION_PLANS_EQUAL( expProvidesPlan, haloExchangePlan.getLocalCommunicationPlan() );
}

/* --------------------------------------------------------------------- */

static HArray<IndexType> randomRequiredIndexes( const Distribution& dist )
{
    const Communicator& comm = dist.getCommunicator();

    const IndexType K = 5 + comm.getRank();   // number of required indexes, different on each proc

    HArray<IndexType> requiredIndexes( K );
    utilskernel::HArrayUtils::setRandom( requiredIndexes, dist.getGlobalSize() - 1 );

    return requiredIndexes;
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildTest )
{
    const IndexType N = 100;     // size of the global array
    
    auto dist = blockDistribution( N );  // block distribution on all available procs

    const Communicator& comm = dist->getCommunicator();

    auto requiredIndexes = randomRequiredIndexes( *dist );

    // Now build a halo exchange plan

    auto plan = haloExchangePlanByRequiredIndexes( requiredIndexes, *dist );

    SCAI_LOG_INFO( logger, comm << ": halo exchange plan = " << plan )

    // Test local indexes for valid entires

    for ( const IndexType localIndex : hostReadAccess( plan.getLocalIndexes() ) )
    {
        BOOST_CHECK( localIndex < dist->getLocalSize() );
    }

    HArray<IndexType> haloIndexes;

    plan.global2HaloV( haloIndexes, requiredIndexes );

    // Make sure that all required indexes are in the halo

    for ( const IndexType haloIndex : hostReadAccess( haloIndexes ) )
    {
        BOOST_REQUIRE( haloIndex != invalidIndex );
    }

    // make sure that translation back works correctly 

    plan.halo2GlobalV( haloIndexes, haloIndexes );
   
    BOOST_TEST( hostReadAccess( haloIndexes ) == hostReadAccess( requiredIndexes ), per_element() );

    // check for valid communication plans

    const CommunicationPlan& localCommPlan = plan.getLocalCommunicationPlan();

    CommunicationPlan expLocalCommPlan = comm.transpose( plan.getHaloCommunicationPlan() );

    CHECK_COMMUNICATION_PLANS_EQUAL( localCommPlan, expLocalCommPlan )
}

/* --------------------------------------------------------------------- */

/** Simple function to set a value in a distributed array. 
 *  It is used to verify correct values in the halo.
 */
template<typename ValueType>
static ValueType globalValue( IndexType globalIndex )
{
    return static_cast<ValueType>( 2 * globalIndex + 1 );
}

template<typename ValueType>
static HArray<ValueType> distributedArray( const Distribution& dist )
{
    HArray<ValueType> localArray;  // local part of the distributed 'global' array

    // use own scope for write access to make sure that access is closed before return

    {   
        IndexType localIndex = 0;   // running local index

        for ( auto& entry : hostWriteOnlyAccess( localArray, dist.getLocalSize() ) )
        {
            entry = globalValue<ValueType>( dist.local2Global( localIndex++ ) );
        }

    }  // filled the local array with 'global' values

    return localArray;    // each processor gets its local part
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( updateTest )
{
    typedef DefaultReal ValueType;

    const IndexType N = 100;     // size of the global array
    
    auto dist = blockDistribution( N );  // block distribution on all available procs

    const Communicator& comm = dist->getCommunicator();

    auto requiredIndexes = randomRequiredIndexes( *dist );

    auto localArray = distributedArray<ValueType>( *dist );

    // Now build a halo exchange plan

    auto plan = haloExchangePlanByRequiredIndexes( requiredIndexes, *dist );

    HArray<ValueType> haloArray;
 
    SCAI_LOG_DEBUG( logger, comm << ": now update halo " )

    plan.updateHalo( haloArray, localArray, comm );

    BOOST_CHECK_EQUAL( haloArray.size(), plan.getHaloSize() );

    SCAI_LOG_INFO( logger, "check updated halo, size = " << plan.getHaloSize() )

    // Now verify that we have received all required values correctly in the halo array

    {
        auto rHalo = hostReadAccess( haloArray );

        for ( const auto& globalIndex : hostReadAccess( requiredIndexes ) )
        {
            IndexType haloIndex = plan.global2Halo( globalIndex );
            BOOST_REQUIRE( haloIndex != invalidIndex );

            ValueType val = rHalo[ haloIndex ];
            ValueType expected = globalValue<ValueType>( globalIndex );
            BOOST_CHECK_EQUAL( val, expected );
        }
    }

    HArray<ValueType> haloArray1;
    HArray<ValueType> tmpArray;        // use this array for sending values
    plan.updateHalo( haloArray1, localArray, comm, tmpArray );

    BOOST_TEST( hostReadAccess( haloArray ) == hostReadAccess( haloArray1 ), per_element() );

    HArray<ValueType> haloArray2;
    plan.updateHaloDirect( haloArray2, tmpArray, comm );
    BOOST_TEST( hostReadAccess( haloArray ) == hostReadAccess( haloArray2 ), per_element() );

    HArray<ValueType> haloArray3;
    std::unique_ptr<tasking::SyncToken> token1(
       plan.updateHaloAsync( haloArray3, localArray, comm ) );
 
    HArray<ValueType> haloArray4;
    std::unique_ptr<tasking::SyncToken> token2(
       plan.updateHaloDirectAsync( haloArray4, tmpArray, comm ) );

    token1->wait();

    BOOST_TEST( hostReadAccess( haloArray ) == hostReadAccess( haloArray3 ), per_element() );

    token2->wait();

    BOOST_TEST( hostReadAccess( haloArray ) == hostReadAccess( haloArray4 ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( updateByTest )
{
    typedef DefaultReal ValueType;

    const IndexType N = 100;     // size of the global array
    
    auto dist = blockDistribution( N );  // block distribution on all available procs

    const Communicator& comm = dist->getCommunicator();

    HArray<IndexType> requiredIndexes( { 5, 13, N - 8, N - 1 } );

    // Now build a halo exchange plan

    auto plan = haloExchangePlanByRequiredIndexes( requiredIndexes, *dist );

    // set local array to zero

    HArray<ValueType> localArray( dist->getLocalSize(), ValueType( 0 ) );

    HArray<ValueType> haloArray( plan.getHaloSize(), ValueType( 1 ) );
 
    plan.updateByHalo( localArray, haloArray, common::BinaryOp::ADD, comm );

    IndexType count = 0;  // verify that each required index is checked

    {
        auto rLocal = hostReadAccess( localArray );

        ValueType expVal( comm.getSize() );  // every processor has had a halo copy

        for ( const IndexType globalI : hostReadAccess( requiredIndexes ) )
        {
            const IndexType localI = dist->global2Local( globalI );

            if ( localI != invalidIndex )
            {
                BOOST_CHECK_EQUAL( rLocal[localI], expVal );
                count++;
            }
        }
    }

    count = comm.sum( count );
  
    BOOST_CHECK_EQUAL( count, requiredIndexes.size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( moveTest )
{
    const IndexType N = 100;     // size of the global array
    
    auto dist = blockDistribution( N );  // block distribution on all available procs

    auto requiredIndexes = randomRequiredIndexes( *dist );

    auto plan = haloExchangePlanByRequiredIndexes( requiredIndexes, *dist );

    // save the pointer of host incarnation of local indexes 

    const IndexType* localIndexesHostPtr = hostReadAccess( plan.getLocalIndexes() ).get();

    HaloExchangePlan plan1( std::move( plan ) );

    BOOST_CHECK_EQUAL( localIndexesHostPtr, hostReadAccess( plan1.getLocalIndexes() ).get() );

    HaloExchangePlan plan2;
    plan2 = std::move( plan1 );

    BOOST_CHECK_EQUAL( localIndexesHostPtr, hostReadAccess( plan2.getLocalIndexes() ).get() );

    HaloExchangePlan plan3;
    plan3.swap( plan2 );

    BOOST_CHECK_EQUAL( localIndexesHostPtr, hostReadAccess( plan3.getLocalIndexes() ).get() );

    CommunicationPlan haloPlan;
    CommunicationPlan localPlan;
    HArray<IndexType> haloIndexes;
    HArray<IndexType> localIndexes;

    plan3.splitUp( haloIndexes, localIndexes, haloPlan, localPlan );

    BOOST_CHECK_EQUAL( localIndexesHostPtr, hostReadAccess( localIndexes ).get() );

    HaloExchangePlan plan4( haloIndexes, localIndexes, haloPlan, localPlan );
    plan4.clear();
    BOOST_CHECK_EQUAL( plan4.getHaloSize(), 0 );
    BOOST_CHECK_EQUAL( plan4.getLocalIndexes().size(), 0 );

    HaloExchangePlan plan5( haloIndexes, localIndexes, haloPlan, localPlan );
    plan5.purge();
    BOOST_CHECK_EQUAL( plan5.getHaloSize(), 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    auto dist = blockDistribution( 15 );
    auto requiredIndexes = randomRequiredIndexes( *dist );

    auto halo1 = haloExchangePlanByRequiredIndexes( requiredIndexes, *dist );

    HaloExchangePlan halo2( halo1 );  // COPY constructor

    CHECK_COMMUNICATION_PLANS_EQUAL( halo1.getHaloCommunicationPlan(), halo2.getHaloCommunicationPlan() )
    CHECK_COMMUNICATION_PLANS_EQUAL( halo1.getLocalCommunicationPlan(), halo2.getLocalCommunicationPlan() )
    
    BOOST_TEST( hostReadAccess( halo1.getLocalIndexes() ) == hostReadAccess( halo2.getLocalIndexes() ), per_element() );
    BOOST_TEST( hostReadAccess( halo1.getHalo2GlobalIndexes() ) == hostReadAccess( halo2.getHalo2GlobalIndexes() ), per_element() );

    // check that global2Halo via generated map works also fine

    HArray<IndexType> haloIndexes1;
    halo1.global2HaloV( haloIndexes1, requiredIndexes );
    HArray<IndexType> haloIndexes2;
    halo2.global2HaloV( haloIndexes2, requiredIndexes );

    BOOST_TEST( hostReadAccess( haloIndexes1 ) == hostReadAccess( haloIndexes2 ), per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
