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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief Contains the implementation of the class GeneralDistributionTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/utilskernel/LArray.hpp>

using namespace scai;
using namespace dmemo;

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

        // get: rank, rank + size, rank + 2 * size, ...

        utilskernel::HArrayUtils::setSequence( localIndexes, rank, size, elemsPerPartition );

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

    hmemo::HArray<IndexType> localIndexes;

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
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    hmemo::HArray<IndexType> localIndexes;

    IndexType N = 2;

    utilskernel::HArrayUtils::setSequence( localIndexes, rank * N, 1, N );

    GeneralDistribution genDist1( size * N, localIndexes, comm );
    const GeneralDistribution& genDist2 = genDist1;
    GeneralDistribution genDist3( size * N, localIndexes, comm );
    GeneralDistribution genDist4( genDist1 );
    BlockDistribution bdist( size * N, comm );
    GeneralDistribution genDist5( bdist );
    GeneralDistribution genDist6( BlockDistribution( size * ( N + 1 ), comm ) );

    BOOST_CHECK_EQUAL( genDist1, genDist2 );  // pointer equality
    BOOST_CHECK_EQUAL( genDist2, genDist1 );  // pointer equality
    BOOST_CHECK_EQUAL( genDist1, genDist3 );  // same constructor equality
    BOOST_CHECK_EQUAL( genDist3, genDist1 );  // same constructor equality
    BOOST_CHECK_EQUAL( genDist1, genDist4 );  // copy equality
    BOOST_CHECK_EQUAL( genDist4, genDist1 );  // copy equality

    if ( size != 1 )
    {
        BOOST_CHECK( genDist1 != bdist );         // do not compare block dist and general dist
    }

    BOOST_CHECK_EQUAL( genDist1, genDist5 );  // but if block is copied to a general dist
    BOOST_CHECK( genDist1 != genDist6 );      // different global sizes
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
