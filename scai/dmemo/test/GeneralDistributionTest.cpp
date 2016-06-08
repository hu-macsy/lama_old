/**
 * @file GeneralDistributionTest.cpp
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
 * @brief Contains the implementation of the class GeneralDistributionTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>

using namespace scai::dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct GeneralDistributionTestConfig
{
    GeneralDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
        rank = comm->getRank();
        size = comm->getSize();
        elemsPerPartition = 10;
        globalSize = elemsPerPartition * size;

        for ( IndexType k = 0; k < elemsPerPartition; ++k )
        {
            localIndexes.push_back( k * size + rank );
        }

        dist = DistributionPtr( new GeneralDistribution( globalSize, localIndexes, comm ) );
    }

    ~GeneralDistributionTestConfig()
    {
        comm = CommunicatorPtr();
    }

    PartitionId rank;
    PartitionId size;

    IndexType elemsPerPartition;
    IndexType globalSize;

    std::vector<IndexType> localIndexes;

    DistributionPtr dist;
};

BOOST_FIXTURE_TEST_SUITE( GeneralDistributionTest, GeneralDistributionTestConfig )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.GeneralDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( generalSizeTest )
{
    BOOST_CHECK( dist->getLocalSize() == elemsPerPartition );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    std::vector<IndexType> localIndexes;
    DistributionPtr generaldist1( new GeneralDistribution( 1, localIndexes, comm ) );
    DistributionPtr generaldist2( generaldist1 );
    DistributionPtr generaldist3( new GeneralDistribution( 1, localIndexes, comm ) );
    DistributionPtr generaldist4( new GeneralDistribution( 3, localIndexes, comm ) );
    BOOST_CHECK( ( *generaldist1 ).isEqual( *generaldist2 ) );
    BOOST_CHECK( !( *generaldist1 ).isEqual( *generaldist3 ) );
    BOOST_CHECK( !( *generaldist1 ).isEqual( *generaldist4 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
