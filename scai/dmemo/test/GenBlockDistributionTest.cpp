/**
 * @file GenBlockDistributionTest.cpp
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
 * @brief Specific tests for GenBlockDistribution.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/utilskernel.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

static CommunicatorPtr comm;

struct GenBlockDistributionTestConfig
{
    GenBlockDistributionTestConfig()
    {
        comm = Communicator::getCommunicatorPtr();
        size = comm->getSize();
        rank = comm->getRank();
        globalSize = size * ( size + 1 );

        for ( PartitionId p = 0; p < size; ++p )
        {
            IndexType localSize = 2 * ( p + 1 );
            localSizes.push_back( localSize );

            for ( IndexType i = 0; i < localSize; i++ )
            {
                theOwners.push_back( p );
            }
        }

        dist = DistributionPtr( new GenBlockDistribution( globalSize, localSizes, comm ) );
    }

    ~GenBlockDistributionTestConfig()
    {
        comm = CommunicatorPtr();
        dist = DistributionPtr();
    }

    DistributionPtr dist;

    PartitionId size;
    PartitionId rank;
    IndexType globalSize;

    std::vector<IndexType> localSizes;
    std::vector<PartitionId> theOwners;
};

BOOST_FIXTURE_TEST_SUITE( GenBlockDistributionTest, GenBlockDistributionTestConfig )

SCAI_LOG_DEF_LOGGER( logger, "Test.GenBlockDistributionTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( genBlockComputeOwnersTest )
{
    hmemo::HArray<IndexType> indexes;
    hmemo::HArray<PartitionId> owners;

    utilskernel::HArrayUtils::setOrder( indexes, globalSize );

    dist->computeOwners( owners, indexes );

    BOOST_REQUIRE_EQUAL( globalSize, owners.size() );
    BOOST_REQUIRE_EQUAL( globalSize, static_cast<IndexType>( theOwners.size() ) );

    // now check for correct owners

    for ( IndexType i = 0; i < globalSize; i++ )
    {
        BOOST_CHECK_EQUAL( theOwners[i], owners[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( genBlockSizeTest )
{
    // setup vector of sizes
    std::vector<IndexType> mlocalSizes; // block size for each partition

    for ( PartitionId p = 0; p < size; ++p )
    {
        // Partition p gets '(p + 1) * 2' elements
        mlocalSizes.push_back( 2 * ( p + 1 ) );
    }

    IndexType globalSize = size * ( size + 1 );
    GenBlockDistribution dist1( globalSize, mlocalSizes, comm );
    BOOST_CHECK( dist1.getLocalSize() == static_cast<IndexType>( 2 * ( rank + 1 ) ) );
    IndexType lb1, ub1;
    dist1.getLocalRange( lb1, ub1 );
    GenBlockDistribution dist2( globalSize, static_cast<IndexType>( 2 * ( rank + 1 ) ), comm );
    IndexType lb2, ub2;
    dist1.getLocalRange( lb2, ub2 );
    BOOST_CHECK( lb1 == lb2 );
    BOOST_CHECK( ub1 == ub2 );
    SCAI_LOG_INFO( logger, "lb = " << lb1 << ", ub = " << ub1 << ", global = " << globalSize )
    GenBlockDistribution dist3( globalSize, lb1, ub1, true, comm );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    IndexType globalSize = comm->getSize();
    IndexType localSize = 1;

    DistributionPtr genblockdist1( new GenBlockDistribution( globalSize, localSize, comm ) );
    DistributionPtr genblockdist2( genblockdist1 );
    DistributionPtr genblockdist3( new GenBlockDistribution( globalSize, localSize, comm ) );
    DistributionPtr genblockdist4( new GenBlockDistribution( 2 * globalSize, 2 * localSize, comm ) );
    BOOST_CHECK( ( *genblockdist1 ).isEqual( *genblockdist2 ) );
    BOOST_CHECK( ( *genblockdist1 ).isEqual( *genblockdist3 ) );
    BOOST_CHECK( !( *genblockdist1 ).isEqual( *genblockdist4 ) );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
