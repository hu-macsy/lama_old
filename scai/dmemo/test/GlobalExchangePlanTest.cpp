/**
 * @file GlobalExchangePlanTest.cpp
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
 * @brief Test routines for the class GlobalExchangePlan
 * @author Thomas Brandes
 * @date 12.12.2018
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/test/data/test_case.hpp>

#include <scai/dmemo/GlobalExchangePlan.hpp>
#include <scai/dmemo/test/TestMacros.hpp>


using namespace scai;
using namespace dmemo;
using namespace hmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( GlobalExchangePlanTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GlobalExchangePlanTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildTest2 )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // this test runs only for 2 processors

    if ( comm->getSize() != 2 )
    {
        return;
    }

    SCAI_LOG_INFO( logger, *comm << ": buildTest" )

    HArray<PartitionId> targets;
    HArray<IndexType> expPerm;
    CommunicationPlan expSendPlan;
    CommunicationPlan expRecvPlan;

    if ( comm->getRank() == 0 )
    {
        targets = HArray<PartitionId>( { 0, 0, 1, 1, 0, 1 } );
        expPerm = HArray<IndexType>( { 0, 1, 4, 2, 3, 5 } );
        expSendPlan = CommunicationPlan( std::vector<IndexType>( { 3, 3 } ) );
        expRecvPlan = CommunicationPlan( std::vector<IndexType>( { 3, 1 } ) );
    }
    else
    {
        targets = HArray<PartitionId>( { 1, 1, 0, 1 } );
        expPerm = HArray<IndexType>( { 2, 0, 1, 3 } );
        expSendPlan = CommunicationPlan( std::vector<IndexType>( { 1, 3 } ) );
        expRecvPlan = CommunicationPlan( std::vector<IndexType>( { 3, 3 } ) );
    }

    GlobalExchangePlan plan( targets, comm );
   
    HArray<IndexType> perm;
    CommunicationPlan sendPlan;
    CommunicationPlan recvPlan;

    plan.splitUp( perm, sendPlan, recvPlan );
    
    SCAI_LOG_INFO( logger, *comm << ": sendPlan = " << sendPlan << ", recvPlan = " << recvPlan )

    CHECK_COMMUNICATION_PLANS_EQUAL( sendPlan, expSendPlan );
    CHECK_COMMUNICATION_PLANS_EQUAL( recvPlan, expRecvPlan );

    BOOST_TEST( hostReadAccess( perm ) == hostReadAccess( expPerm ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildTest4 )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    // this test runs only for 4 processors

    if ( comm->getSize() != 4 )
    {
        return;
    }

    SCAI_LOG_INFO( logger, *comm << ": buildTest" )

    HArray<PartitionId> targets;
    HArray<IndexType> expPerm;
    CommunicationPlan expSendPlan;
    CommunicationPlan expRecvPlan;

    HArray<IndexType> sendValues;
    HArray<IndexType> expRecvValues;

    // Example from Users Guide for dmemo module

    if ( comm->getRank() == 0 )
    {
        targets = HArray<PartitionId>( { 1, 2, 1, 1, 2 } );
        sendValues = HArray<PartitionId>( { 13, 26, 12, 11, 23 } );

        expPerm = HArray<IndexType>( { 0, 2, 3, 1, 4 } );
        expSendPlan = CommunicationPlan( std::vector<IndexType>( { 0, 3, 2, 0 } ) );
        expRecvPlan = CommunicationPlan( std::vector<IndexType>( { 0, 1, 1, 0 } ) );
        expRecvValues = HArray<PartitionId>( { 5, 8 } );
    }
    else if ( comm->getRank() == 1 )
    {
        targets = HArray<PartitionId>( { 2, 0, 3, 2 } );
        sendValues = HArray<PartitionId>( { 22, 5, 35, 27 } );

        expPerm = HArray<IndexType>( { 1, 0, 3, 2 } );
        expSendPlan = CommunicationPlan( std::vector<IndexType>( { 1, 0, 2, 1 } ) );
        expRecvPlan = CommunicationPlan( std::vector<IndexType>( { 3, 0, 3, 1 } ) );
        expRecvValues = HArray<PartitionId>( { 13, 12, 11, 13, 19, 16, 18  } );
    }
    else if ( comm->getRank() == 2 )
    {
        targets = HArray<PartitionId>( { 3, 0, 1, 3, 1, 1 } );
        sendValues = HArray<PartitionId>( { 30, 8, 13, 31, 19, 16 } );

        expPerm = HArray<IndexType>( { 1, 2, 4, 5, 0, 3  } );
        expSendPlan = CommunicationPlan( std::vector<IndexType>( { 1, 3, 0, 2 } ) );
        expRecvPlan = CommunicationPlan( std::vector<IndexType>( { 2, 2, 0, 2 } ) );
        expRecvValues = HArray<PartitionId>( { 26, 23, 22, 27, 29, 26 } );
    }
    else if ( comm->getRank() == 3 )
    {
        // Attention: here we have added an invalid owner, should only give a warning

        targets = HArray<PartitionId>( { 2, 1, 2, 5 } );
        sendValues = HArray<PartitionId>( { 29, 18, 26, 0 } );

        expPerm = HArray<IndexType>( { 1, 0, 2 } );
        expSendPlan = CommunicationPlan( std::vector<IndexType>( { 0, 1, 2, 0 } ) );
        expRecvPlan = CommunicationPlan( std::vector<IndexType>( { 0, 1, 2, 0 } ) );
        expRecvValues = HArray<PartitionId>( { 35, 30, 31 } );
    }

    GlobalExchangePlan plan( targets, comm );

    HArray<IndexType> recvValues;
    plan.exchange( recvValues, sendValues );

    HArray<IndexType> returnValues( sendValues.size(), IndexType( 0 ) );   // will be the same as sendValues
    plan.exchangeBack( returnValues, recvValues );
   
    HArray<IndexType> perm;
    CommunicationPlan sendPlan;
    CommunicationPlan recvPlan;

    plan.splitUp( perm, sendPlan, recvPlan );
    
    SCAI_LOG_INFO( logger, *comm << ": sendPlan = " << sendPlan << ", recvPlan = " << recvPlan )

    CHECK_COMMUNICATION_PLANS_EQUAL( sendPlan, expSendPlan );
    CHECK_COMMUNICATION_PLANS_EQUAL( recvPlan, expRecvPlan );

    BOOST_TEST( hostReadAccess( perm ) == hostReadAccess( expPerm ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( recvValues ) == hostReadAccess( expRecvValues ), boost::test_tools::per_element() );
    BOOST_TEST( hostReadAccess( returnValues ) == hostReadAccess( sendValues ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
