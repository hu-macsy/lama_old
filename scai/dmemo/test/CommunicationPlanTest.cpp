/**
 * @file CommunicationPlanTest.cpp
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

    CommunicationPlan requiredPlan( reqQuantities.data(), reqQuantities.size(), compressFlag );

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

    CommunicationPlan requiredPlan2( comm->getSize(), &reqOwners[0], static_cast<IndexType>( reqOwners.size() ) );

    BOOST_CHECK_EQUAL( requiredPlan2.totalQuantity(), static_cast<IndexType>( reqOwners.size() ) );

    CommunicationPlan requiredPlan1( &reqQuantities[0], reqQuantities.size() );

    // verify that both plans are same

    BOOST_CHECK_EQUAL( requiredPlan2.totalQuantity(), requiredPlan1.totalQuantity() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeTest )
{
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    std::vector<IndexType> reqQuantities;

    setQuantities( reqQuantities, *comm );

    CommunicationPlan requiredPlan( reqQuantities.data(), reqQuantities.size() );

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

    CommunicationPlan tmpPlan( reqQuantities.data(), reqQuantities.size(), compressed );

    const IndexType nMult = 2;

    CommunicationPlan requiredPlan( tmpPlan, nMult );

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

    CommunicationPlan requiredPlan( reqQuantities.data(), reqQuantities.size() );

    CommunicationPlan providesPlan;

    providesPlan.allocateTranspose( requiredPlan, *comm );

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

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
