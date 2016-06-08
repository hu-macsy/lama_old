/**
 * @file MetisDistributionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
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
 * @endlicense
 *
 * @brief MetisDistributionTest.cpp
 * @author Lauretta Schubert
 * @date 01.07.2013
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/MetisDistribution.hpp>
#include <scai/dmemo.hpp>

#include <scai/common/shared_ptr.hpp>

#include <vector>

using namespace scai::dmemo;

using scai::common::shared_ptr;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct MetisDistributionTestConfig
{
    MetisDistributionTestConfig()
    {
        comm = Communicator::getCommunicator();

        rank = comm->getRank();
        size = comm->getSize();

        globalSize = 17;

        DistributionPtr d ( new NoDistribution( globalSize ) );

        matrix.reset( new Distributed( d ) );

        // weights

        float weight = static_cast<float>( 1.0 / size );
        parts.reserve( size );

        for ( int i = 0; i < size - 1; ++i )
        {
            parts[i] = weight;
        }

        parts[ size - 1 ] = 1.0f - ( size - 1 ) * weight;
        dist = DistributionPtr( new MetisDistribution( comm, *matrix, parts ) );
    }

    ~MetisDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType globalSize;

    DistributionPtr dist;

    shared_ptr<Distributed> matrix;

    std::vector<float> parts;
};

BOOST_FIXTURE_TEST_SUITE( MetisDistributionTest, MetisDistributionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.MetisDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    Distributed distMatrix( *matrix );
    DistributionPtr generaldist1( new MetisDistribution( comm, distMatrix, parts ) );
    DistributionPtr generaldist2( generaldist1 );
    DistributionPtr generaldist3( new MetisDistribution( comm, distMatrix, parts ) );
    BOOST_CHECK(  ( *generaldist1 ).isEqual( *generaldist2 ) );
    BOOST_CHECK( !( *generaldist1 ).isEqual( *generaldist3 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
