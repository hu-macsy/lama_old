/**
 * @file CommunicationPlanTest.cpp
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
 * @brief Test methods for the class CommunicationPlan
 * @author Thomas Brandes
 * @date 09.05.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/test/TestMacros.hpp>

#include <scai/logging.hpp>

#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/CommunicationPlan.hpp>

using namespace scai;
using namespace dmemo;

using hmemo::HArray;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommunicationPlanTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CommunicationPlanTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( generalTest )
{
    std::vector<IndexType> quantities( { 0, 0, 2, 3, 0, 1 } );

    CommunicationPlan plan( quantities );

    BOOST_REQUIRE_EQUAL( plan.size(), 3 );   // there are exactly 3 non-zero entries

    BOOST_CHECK_EQUAL( plan.maxQuantity(), 3 );
    BOOST_CHECK_EQUAL( plan.totalQuantity(), 6 );

    // check the second entry, partner = 3, offset = 2

    const CommunicationPlan::Entry& entry = plan[1];

    BOOST_CHECK_EQUAL( entry.partitionId, PartitionId( 3 ) );
    BOOST_CHECK_EQUAL( entry.offset, 2 );

    const CommunicationPlan::Entry& entry2 = plan[2];

    BOOST_CHECK_EQUAL( entry2.partitionId, PartitionId( 5 ) );
    BOOST_CHECK_EQUAL( entry2.offset, 5 );

    IndexType size;
    IndexType offset;

    plan.getInfo( size, offset, 1 );
    BOOST_CHECK_EQUAL( size, IndexType( 0 ) );
    plan.getInfo( size, offset, 3 );
    BOOST_CHECK_EQUAL( size, IndexType( 3 ) );
    BOOST_CHECK_EQUAL( offset, IndexType( 2 ) );
}

/* --------------------------------------------------------------------- */

// Provide for test a closed formula that delivers required quantities

static IndexType required( const PartitionId rankRequires, const PartitionId rankProvides )
{
    if ( rankRequires == rankProvides )
    {
        return 0;
    }

    return 2 * rankProvides + rankRequires;
}

static std::vector<IndexType> getRequiredSizes( const Communicator& comm )
{
    std::vector<IndexType> quantities;

    PartitionId rank = comm.getRank();
    PartitionId size = comm.getSize();

    quantities.reserve( size );

    for ( PartitionId p = 0; p < size; ++p )
    {
        // set number of entries this processor ( rank ) requires from other processor p

        quantities.push_back( required( rank, p ) );
    }

    return quantities;
}

static std::vector<IndexType> getProvidesSizes( const Communicator& comm )
{
    std::vector<IndexType> quantities;

    PartitionId rank = comm.getRank();
    PartitionId size = comm.getSize();

    quantities.resize( size );

    for ( PartitionId p = 0; p < size; ++p )
    {
        // set number of entries this processor ( rank ) provides to other processor p

        quantities[p] = required( p, rank );
    }

    return quantities;
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocatePlanTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    CommunicationPlan requiredPlan( getRequiredSizes( *comm ) );

    // verify that requiredPlan is correctly set up

    IndexType offsetCheck = 0;

    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        IndexType n = requiredPlan[p].quantity;
        PartitionId partitionId = requiredPlan[p].partitionId;
        IndexType nExpected = required( rank, partitionId );
        BOOST_CHECK_EQUAL( n, nExpected );
        BOOST_CHECK_EQUAL( requiredPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( offsetCheck, requiredPlan.totalQuantity() );

    requiredPlan.clear();

    BOOST_CHECK_EQUAL( IndexType( 0 ), requiredPlan.size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    CommunicationPlan requiredPlan( getRequiredSizes( *comm ) );

    std::ostringstream out;

    out << requiredPlan;

    BOOST_CHECK( out.str().length() > 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    CommunicationPlan plan1( getRequiredSizes( *comm ) );
    CommunicationPlan plan2( plan1 );

    CHECK_COMMUNICATION_PLANS_EQUAL( plan1, plan2 )

    CommunicationPlan plan3;   // zero plan
    plan3 = plan1;

    CHECK_COMMUNICATION_PLANS_EQUAL( plan1, plan3 )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( moveTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<IndexType> sizes( { 3, 4, 1, 0 } );
    CommunicationPlan plan1( sizes );

    // we test for correct move by checking for pointer in the allocated data structures

    const CommunicationPlan::Entry& entry1 = plan1[1];

    CommunicationPlan plan2( std::move( plan1 ) );
    const CommunicationPlan::Entry& entry2 = plan2[1];

    BOOST_CHECK_EQUAL( &entry1, &entry2 );

    CommunicationPlan plan3 = std::move( plan2 );
    const CommunicationPlan::Entry& entry3 = plan3[1];

    BOOST_CHECK_EQUAL( &entry1, &entry3 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( multiplyTest )
{
    const IndexType nMult = 2;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<IndexType> sizes = getRequiredSizes( *comm );

    CommunicationPlan plan1 ( sizes );
    plan1.multiplyConst( nMult );

    // compare it against the alternative by building it with new sizes

    for ( auto& entry : sizes ) 
    {
        entry *= nMult;
    }

    CommunicationPlan plan2( sizes );

    CHECK_COMMUNICATION_PLANS_EQUAL( plan1, plan2 )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructNTest )
{
    const IndexType N = 3;

    std::vector<IndexType> commSizes( { 0, 2, 1, 0, 3 } );
    std::vector<IndexType> commSizesN( { 0 * N, 2 * N, 1 * N, 0 * N, 3 * N } );

    CommunicationPlan plan( commSizes );
    CommunicationPlan planExpected( commSizesN );

    auto planN = plan.constructN( N );

    CHECK_COMMUNICATION_PLANS_EQUAL( planExpected, planN )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( raggedTest )
{
    HArray<IndexType> commSizes( { 0, 2, 1, 0, 3 } );      // sizes for original plan

    HArray<IndexType> raggedSizes( { 2, 1, 3, 2, 0, 3 } );  // array of array sizes
    HArray<IndexType> raggedOffsets( { 0, 2, 3, 6, 8, 8, 11 } );  // = sizes2Offsets( raggedSizes )

    HArray<IndexType> raggedCommSizes( { 0, 2 + 1, 3, 0, 2 + 0 + 3 } );

    CommunicationPlan plan( commSizes );
    CommunicationPlan planExpected( raggedCommSizes );

    auto plan1 = plan.constructRaggedBySizes( raggedSizes );
    auto plan2 = plan.constructRaggedByOffsets( raggedOffsets );

    CHECK_COMMUNICATION_PLANS_EQUAL( planExpected, plan1 )
    CHECK_COMMUNICATION_PLANS_EQUAL( planExpected, plan2 )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTransposeTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    CommunicationPlan requiredPlan( getRequiredSizes( *comm ) );

    auto providesPlan = comm->transpose( requiredPlan );

    CommunicationPlan expectedProvidesPlan( getProvidesSizes( *comm ) );

    CHECK_COMMUNICATION_PLANS_EQUAL( providesPlan, expectedProvidesPlan )
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( singleEntryTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType quantities[] = { 3, 0, 5, 2 };

    CommunicationPlan plan( quantities, 4 );
   
    BOOST_CHECK_EQUAL( plan.size(), 3 );
    BOOST_CHECK_EQUAL( plan.maxQuantity(), 5 );

    plan.defineBySingleEntry( 4, 2 );

    BOOST_CHECK_EQUAL( plan.size(), 1 );

    BOOST_CHECK_EQUAL( plan.maxQuantity(), 4 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getInfoTest )
{
    // Idea: build a communication plan and call getInfo where quantity and offset are known

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    CommunicationPlan requiredPlan( getRequiredSizes( *comm ) );

    IndexType expectedOffset = 0;

    for ( PartitionId p = 0; p < comm->getSize(); ++p )
    {
        IndexType n;
        IndexType offset;
        requiredPlan.getInfo( n, offset, p );

        BOOST_CHECK_EQUAL( n, required( rank, p ) );

        if ( n > 0 )
        {
            BOOST_CHECK_EQUAL( expectedOffset, offset );
            expectedOffset += n;
        }
    }

    BOOST_CHECK_EQUAL( expectedOffset, requiredPlan.totalQuantity() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( extractPlanTest )
{
    // Idea: build a communication plan and call getInfo where quantity and offset are known

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    CommunicationPlan requiredPlan( getRequiredSizes( *comm ) );

    CommunicationPlan singlePlan;

    for ( PartitionId p = 0; p < comm->getSize(); ++p )
    {
        singlePlan.extractPlan( requiredPlan, p );

        if ( singlePlan.size() == 0 )
        {
            BOOST_CHECK_EQUAL( IndexType( 0 ), required( rank, p ) );
        }
        else
        {
            // singlePlan must have exactly one entry for this p

            BOOST_CHECK_EQUAL( singlePlan.size(), IndexType( 1 ) );
            BOOST_CHECK_EQUAL( p, singlePlan[0].partitionId );
            BOOST_CHECK_EQUAL( required( rank, p ), singlePlan[0].quantity );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
