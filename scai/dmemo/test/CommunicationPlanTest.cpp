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

#include <scai/logging.hpp>

#include <scai/dmemo/Communicator.hpp>
#include <scai/dmemo/CommunicationPlan.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommunicationPlanTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CommunicationPlanTest" )

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

static void setQuantities( std::vector<IndexType>& quantities, const Communicator& comm )
{
    PartitionId rank = comm.getRank();
    PartitionId size = comm.getSize();

    quantities.resize( size );

    for ( PartitionId p = 0; p < size; ++p )
    {
        // set number of entries this processor ( rank ) requires from other processor p

        quantities[p] = required( rank, p );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocatePlanTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    bool compressFlag = false;

    auto requiredPlan = CommunicationPlan::buildBySizes( reqQuantities.data(), reqQuantities.size(), compressFlag );

    BOOST_CHECK( !requiredPlan.compressed() );

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

BOOST_AUTO_TEST_CASE( allocatePlanByOffsetTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    std::vector<IndexType> quantities;

    setQuantities( quantities, *comm );

    std::vector<IndexType> offsets;

    IndexType offs = 0;

    offsets.push_back( offs );

    for ( size_t i = 0; i < quantities.size(); ++i )
    {
        offs += quantities[i];
        offsets.push_back( offs );
    }

    BOOST_REQUIRE_EQUAL( offsets.size(), quantities.size() + 1 );

    bool compressFlag = true;

    CommunicationPlan requiredPlan;

    requiredPlan.allocateByOffsets( &offsets[0], quantities.size(), compressFlag );

    BOOST_CHECK( requiredPlan.compressed() );

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

    BOOST_CHECK_EQUAL( offs, offsetCheck );
    BOOST_CHECK_EQUAL( offsetCheck, requiredPlan.totalQuantity() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    std::vector<PartitionId> reqOwners;

    for ( PartitionId owner = 0; owner < static_cast<PartitionId>( reqQuantities.size() ); ++owner )
    {
        for ( IndexType k = 0; k < reqQuantities[owner]; ++k )
        {
            reqOwners.push_back( owner );
        }
    }

    const IndexType* reqQuantitiesBegin = reqQuantities.size() > 0 ? &reqQuantities[0] : NULL;

    const IndexType* reqOwnersBegin = reqOwners.size() > 0 ? &reqOwners[0] : NULL;

    CommunicationPlan requiredPlan2( comm->getSize(), reqOwnersBegin, static_cast<IndexType>( reqOwners.size() ) );

    BOOST_CHECK_EQUAL( requiredPlan2.totalQuantity(), static_cast<IndexType>( reqOwners.size() ) );

    auto requiredPlan1 = CommunicationPlan::buildBySizes( reqQuantitiesBegin, reqQuantities.size() );

    // verify that both plans are same

    BOOST_CHECK_EQUAL( requiredPlan2.totalQuantity(), requiredPlan1.totalQuantity() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    auto requiredPlan = CommunicationPlan::buildBySizes( reqQuantities.data(), reqQuantities.size() );

    std::ostringstream out;

    out << requiredPlan;

    BOOST_CHECK( out.str().length() > 0 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    bool compressed = true;

    auto tmpPlan = CommunicationPlan::buildBySizes( reqQuantities.data(), reqQuantities.size(), compressed );

    const IndexType nMult = 2;

    CommunicationPlan requiredPlan( tmpPlan );
    requiredPlan.multiplyConst( nMult );

    BOOST_CHECK( requiredPlan.allocated() );
    BOOST_CHECK( requiredPlan.compressed() );

    // verify that requiredPlan is correctly set up

    IndexType offsetCheck = 0;

    for ( PartitionId p = 0; p < requiredPlan.size(); ++p )
    {
        IndexType n = requiredPlan[p].quantity;
        PartitionId partitionId = requiredPlan[p].partitionId;
        IndexType nExpected = required( rank, partitionId );
        BOOST_CHECK_EQUAL( n, nMult * nExpected );
        BOOST_CHECK_EQUAL( requiredPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( offsetCheck, requiredPlan.totalQuantity() );

    requiredPlan.purge();

    BOOST_CHECK_EQUAL( IndexType( 0 ), requiredPlan.size() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( allocateTransposeTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    auto requiredPlan = CommunicationPlan::buildBySizes( reqQuantities.data(), reqQuantities.size() );
    auto providesPlan = requiredPlan.transpose( *comm );

    IndexType offsetCheck = 0;

    for ( PartitionId p = 0; p < providesPlan.size(); ++p )
    {
        IndexType n = providesPlan[p].quantity;
        PartitionId partitionId = providesPlan[p].partitionId;
        IndexType nExpected = required( partitionId, rank );   // here switched arguments
        BOOST_CHECK_EQUAL( n, nExpected );
        BOOST_CHECK_EQUAL( providesPlan[p].offset, offsetCheck );
        offsetCheck += n;
    }

    BOOST_CHECK_EQUAL( offsetCheck, providesPlan.totalQuantity() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( singleEntryTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    const IndexType quantities[] = { 3, 0, 5, 2 };

    auto plan = CommunicationPlan::buildBySizes( quantities, 4, true );
   
    BOOST_CHECK_EQUAL( plan.size(), 3 );
    BOOST_CHECK_EQUAL( plan.maxQuantity(), 5 );

    plan.singleEntry( 2, 4 );

    BOOST_CHECK_EQUAL( plan.size(), 1 );

    BOOST_CHECK_EQUAL( plan.maxQuantity(), 4 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getInfoTest )
{
    // Idea: build a communication plan and call getInfo where quantity and offset are known

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    bool compressFlag = true;  // make sure that p in getInfo is not same as entry pos

    auto requiredPlan = CommunicationPlan::buildBySizes( reqQuantities.data(), reqQuantities.size(), compressFlag );

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

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    bool compressFlag = true;  // make sure that p in getInfo is not same as entry pos

    auto requiredPlan = CommunicationPlan::buildBySizes( reqQuantities.data(), reqQuantities.size(), compressFlag );

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
