/**
 * @file NoDistributionTest.cpp
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
 * @brief Specific tests for the class NoDistributionTest derived from Distribution.
 * @author Thmas Brandes
 * @date 01.08.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo/NoDistribution.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct NoDistributionTestConfig
{
    NoDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
        rank = comm->getRank();
        size = comm->getSize();
        blockSize = 17;
        dist = DistributionPtr( new NoDistribution( blockSize * size ) );
    }

    ~NoDistributionTestConfig()
    {
        comm = CommunicatorPtr();
        dist = DistributionPtr();
    }

    PartitionId rank;
    PartitionId size;
    IndexType blockSize;

    std::vector<IndexType> nonLocalIndexes;

    DistributionPtr dist;
};

BOOST_FIXTURE_TEST_SUITE( NoDistributionTest, NoDistributionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.NoDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( CtorTest )
{
    const IndexType n = 20;

    DistributionPtr nodist( new NoDistribution( n ) );
    BOOST_CHECK_EQUAL( nodist->getCommunicatorPtr()->getType(), CommunicatorType::NO );
    BOOST_CHECK_EQUAL( nodist->getGlobalSize(), n );
    BOOST_CHECK_EQUAL( nodist->getLocalSize(), n );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isLocalTest )
{
    BOOST_CHECK( dist->isLocal( 0 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    DistributionPtr nodist1( new NoDistribution( 1 ) );
    DistributionPtr nodist2( nodist1 );
    DistributionPtr nodist3( new NoDistribution( 1 ) );
    DistributionPtr nodist4( new NoDistribution( 2 ) );
    BOOST_CHECK( ( *nodist1 ).isEqual( *nodist2 ) );
    BOOST_CHECK( ( *nodist1 ).isEqual( *nodist3 ) );
    BOOST_CHECK( !( *nodist1 ).isEqual( *nodist4 ) );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
