/**
 * @file GenBlockDistributionTest.cpp
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

        dist = genBlockDistributionBySize( 2 * ( rank + 1 ), comm );
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

    utilskernel::HArrayUtils::setOrder( indexes, globalSize );

    auto owners = dist->owner( indexes );

    BOOST_TEST( hostReadAccess( owners ) == theOwners, boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( genBlockByWeightTest )
{
    const IndexType N = 24;
  
    auto comm = Communicator::getCommunicatorPtr();

    bool  even   = comm->getRank() % 2 == 0;
    float weight = even ? 1.0f : 2.0f ;

    if ( comm->getRank() > 3 ) 
    {
        weight = 0.0f;
    }

    auto dist = genBlockDistributionByWeight( N, weight, comm );

    if ( comm->getSize() == 1 )
    {
        IndexType expLocalSize = N;
        BOOST_CHECK_EQUAL( expLocalSize, dist->getLocalSize() );
    }
    else if ( comm->getSize() == 2 )
    {
        IndexType expLocalSize = even ? N / 3 : 2 * N / 3; 
        BOOST_CHECK_EQUAL( expLocalSize, dist->getLocalSize() );
    }
    else if ( comm->getSize() == 3 )
    {
        IndexType expLocalSize = even ?  N / 4 : 2 * N / 4;
        BOOST_CHECK_EQUAL( expLocalSize, dist->getLocalSize() );
    }
    else 
    {
        IndexType expLocalSize = even ?  N / 6 : 2 * N / 6;

        if ( comm->getRank() > 3 )
        {
            expLocalSize = 0;
        }

        SCAI_LOG_DEBUG( logger, *comm << ": exp local size = " << expLocalSize << " for dist = " << *dist )

        BOOST_CHECK_EQUAL( expLocalSize, dist->getLocalSize() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructByOffsetTest )
{
    const IndexType N = comm->getSize();
    const IndexType R = comm->getRank();

    const IndexType offset = R;
    auto dist = genBlockDistributionByOffset( N, offset, comm );
   
    BOOST_CHECK_EQUAL( dist->getLocalSize(), IndexType( 1 ) );
    BOOST_CHECK_EQUAL( dist->lb(), offset );

    const IndexType offset1 = R % 2 == 0 ? invalidIndex : R;
    auto dist1 = genBlockDistributionByOffset( N, offset1, comm );
 
    BOOST_CHECK_EQUAL( comm->sum( dist1->getLocalSize() ), N );

    BOOST_CHECK_THROW(
    {
        auto dist2 = genBlockDistributionByOffset( N, N + 1, comm );
    }, common::Exception );
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

    auto dist1 = genBlockDistributionBySizes( mlocalSizes, comm );
    auto dist2 = genBlockDistributionBySize( static_cast<IndexType>( 2 * ( rank + 1 ) ), comm );

    BOOST_CHECK_EQUAL( dist1->getGlobalSize(), globalSize );
    BOOST_CHECK_EQUAL( dist2->getGlobalSize(), globalSize );
    BOOST_CHECK_EQUAL( dist1->getLocalSize(), static_cast<IndexType>( 2 * ( rank + 1 ) ) );
    BOOST_CHECK_EQUAL( dist1->getLocalSize(), dist2->getLocalSize() );

    IndexType lb1 = dist1->lb();
    IndexType ub1 = dist1->ub();
    IndexType lb2 = dist2->lb();
    IndexType ub2 = dist2->ub();
    BOOST_CHECK_EQUAL( lb1, lb2 );
    BOOST_CHECK_EQUAL( ub1, ub2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    IndexType localSize = 1;

    auto genblockdist1 = genBlockDistributionBySize( localSize, comm );
    DistributionPtr genblockdist2( genblockdist1 );
    auto genblockdist3 = genBlockDistributionBySize( localSize, comm );
    auto genblockdist4 = genBlockDistributionBySize( 2 * localSize, comm );
    BOOST_CHECK( ( *genblockdist1 ).isEqual( *genblockdist2 ) );
    BOOST_CHECK( ( *genblockdist1 ).isEqual( *genblockdist3 ) );
    BOOST_CHECK( !( *genblockdist1 ).isEqual( *genblockdist4 ) );

    BOOST_CHECK( genblockdist1->isSameGenBlockDistribution( *genblockdist3 ) );
    BOOST_CHECK( genblockdist1->isBlockDistribution() );
    BOOST_CHECK( genblockdist4->isBlockDistribution() );

    const IndexType N = 1000 * comm->getSize();
    auto genblockdist5 = genBlockDistributionByWeight( N, 1.0f );
    BOOST_CHECK( genblockdist5->isBlockDistribution() );
    auto genblockdist6 = genBlockDistributionByWeight( N, float( comm->getRank() + 1 ) );
    if ( comm->getSize() > 1 )
    {
        BOOST_CHECK( !genblockdist6->isBlockDistribution() );
    }
    auto genblockdist7 = genBlockDistributionByWeight( N, 2.0f );
    BOOST_CHECK( genblockdist7->isSameGenBlockDistribution( *genblockdist5 ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
