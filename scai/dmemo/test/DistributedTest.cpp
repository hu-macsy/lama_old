/**
 * @file DistributedTest.cpp
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
 * @brief Test for base class Distributed via Mock class
 * @author Thomas Brandes
 * @date 21.07.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/Distributed.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>
#include <scai/utilskernel.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DistributedTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DistributedTest" )

/* --------------------------------------------------------------------- */

class MockDistributed : public Distributed
{
public:

    MockDistributed( DistributionPtr dist ) : Distributed( dist )
    {
    }

    void redistribute( DistributionPtr newDist )
    {
        setDistributionPtr( newDist );
    }

    void swap( MockDistributed& other )
    {
        Distributed::swap( other );
    }
};

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    const IndexType globalSize = 5;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr dist( new BlockDistribution( globalSize, comm ) );

    MockDistributed mock( dist );

    BOOST_CHECK_EQUAL( mock.getDistributionPtr(), dist );
    BOOST_CHECK_EQUAL( mock.getDistribution(), *dist );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( copyTest )
{
    const IndexType globalSize = 5;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr dist( new CyclicDistribution( globalSize, 1, comm ) );

    MockDistributed mock1( dist );

    // default copy constructor of MockDistributedcalls copy constructor of Distributed

    MockDistributed mock2( mock1 );

    BOOST_CHECK_EQUAL( mock1.getDistribution(), mock2.getDistribution() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( swapTest )
{
    const IndexType globalSize = 5;

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    DistributionPtr dist1( new BlockDistribution( globalSize, comm ) );
    DistributionPtr dist2( new CyclicDistribution( globalSize, 1, comm ) );

    MockDistributed mock1( dist1 );
    MockDistributed mock2( dist1 );
    mock2.redistribute( dist2 );

    BOOST_CHECK_EQUAL( mock1.getDistribution(), *dist1 );
    BOOST_CHECK_EQUAL( mock2.getDistribution(), *dist2 );

    mock1.swap( mock2 );

    BOOST_CHECK_EQUAL( mock1.getDistribution(), *dist2 );
    BOOST_CHECK_EQUAL( mock2.getDistribution(), *dist1 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */

