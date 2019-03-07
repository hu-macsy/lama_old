/**
 * @file GeneralDistributionTest.cpp
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
 * @brief Specific tests for derived distribution class GeneralDistributionTest.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/dmemo/BlockDistribution.hpp>
#include <scai/dmemo/CyclicDistribution.hpp>

#include <scai/utilskernel.hpp>

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

        const IndexType first = static_cast<IndexType>( rank );
        const IndexType inc   = static_cast<IndexType>( size );

        utilskernel::HArrayUtils::setSequence( localIndexes, first, inc, elemsPerPartition );

        dist = DistributionPtr( new GeneralDistribution( globalSize, localIndexes, false, comm ) );
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

std::shared_ptr<const GeneralDistribution> buildCyclic( 
    const IndexType elemsPerProcessor, 
    CommunicatorPtr comm = Communicator::getCommunicatorPtr() )

{
    const IndexType rank = comm->getRank();
    const IndexType NP   = comm->getSize();

    hmemo::HArray<IndexType> localIndexes;
    utilskernel::HArrayUtils::setSequence( localIndexes, rank, NP, elemsPerProcessor );

    return generalDistributionUnchecked( elemsPerProcessor * NP, std::move( localIndexes ), comm );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( generalSizeTest )
{
    const IndexType elemsPerProcessor = 10;

    auto dist = buildCyclic( elemsPerProcessor );

    const Communicator& comm = dist->getCommunicator();

    BOOST_CHECK_EQUAL( dist->getLocalSize(), elemsPerProcessor );
    BOOST_CHECK_EQUAL( dist->getGlobalSize(), comm.sum( dist->getLocalSize() ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( isEqualTest )
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    hmemo::HArray<IndexType> localIndexes;

    IndexType N = 2;

    const IndexType first = rank * N;
    const IndexType inc   = 1;

    bool checkFlag = true;   // not really need but for convenience

    utilskernel::HArrayUtils::setSequence( localIndexes, first, inc, N );

    GeneralDistribution genDist1( size * N, localIndexes, checkFlag, comm );
    const GeneralDistribution& genDist2 = genDist1;
    GeneralDistribution genDist3( size * N, localIndexes, checkFlag, comm );

    // general distributions can be compared with each other

    BOOST_CHECK_EQUAL( genDist1, genDist2 );  // pointer equality
    BOOST_CHECK_EQUAL( genDist2, genDist1 );  // pointer equality
    BOOST_CHECK_EQUAL( genDist1, genDist3 );  // same constructor equality
    BOOST_CHECK_EQUAL( genDist3, genDist1 );  // same constructor equality
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getBlockDistributedOwnersTest )
{
    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    const GeneralDistribution* genDist = dynamic_cast<const GeneralDistribution*>( dist.get() );

    BOOST_REQUIRE( genDist != NULL );

    const hmemo::HArray<PartitionId>& localOwners = genDist->getMyBlockDistributedOwners();

    hmemo::HArray<PartitionId> expectedOwners;

    {
        auto wExpected = hmemo::hostWriteOnlyAccess( expectedOwners, elemsPerPartition );

        for ( IndexType i = 0; i < elemsPerPartition; ++i )
        {
            IndexType globalIndex = i + rank * elemsPerPartition;   
            PartitionId owner = globalIndex % size; 
            wExpected[i] = owner;
        }
    }
  
    BOOST_TEST( hostReadAccess( expectedOwners ) == hostReadAccess( localOwners ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( redistConstructorTest )
{
    IndexType N = 15;

    TestDistributions allDist( N );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        const Distribution& dist = *allDist[i];

        if ( dist.getCommunicator() != *comm )
        {
            // dist is NoDistribution, cannot be redistributed

            continue;
        }

        IndexType nLocal = dist.getLocalSize();

        hmemo::HArray<PartitionId> owners;  // new owners for my local indexes

        {
            hmemo::WriteOnlyAccess<PartitionId> wOwners( owners, nLocal );

            for ( IndexType i = 0; i < nLocal; ++i )
            {
                // choose owner as if it will be a Cyclic(1) distribution

                IndexType globalIndex = dist.local2Global( i );
                wOwners[i] = globalIndex % size;
            }
        }

        SCAI_LOG_DEBUG( logger, "redistribute, dist = " << dist << ", owners = " << owners )

        auto gdist = generalDistributionByNewOwners( dist, owners );
        auto cdist = cyclicDistribution( N, 1, comm );

        BOOST_REQUIRE_EQUAL( gdist->getLocalSize(), cdist->getLocalSize() );

        for ( IndexType i = 0; i < cdist->getLocalSize(); ++i )
        {
            BOOST_CHECK_EQUAL( cdist->local2Global( i ), gdist->local2Global( i ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildByOwnersTest )
{
    auto comm = Communicator::getCommunicatorPtr();

    const IndexType N = 15;

    auto size = comm->getSize();
    auto rank = comm->getRank();

    hmemo::HArray<PartitionId> owners;

    decltype( rank ) root = 0;

    if ( rank == root )
    {
        // define cylcic(1) ownership

        auto wOwners = hostWriteOnlyAccess( owners, N );  

        for ( IndexType i = 0; i < N; ++i )
        {
            wOwners[i] = i % size;
        }
    }

    auto dist = generalDistributionBySingleOwners( owners, root, comm );

    // now prove that each processor has the right local values

    hmemo::HArray<IndexType> myIndexes;

    dist->getOwnedIndexes( myIndexes );

    for ( auto myIndex : hostReadAccess( myIndexes ) )
    {
        BOOST_CHECK_EQUAL( static_cast<decltype( rank )>( myIndex % size ), rank );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildByOwnersFailTest )
{
    auto comm = Communicator::getCommunicatorPtr();

    auto size = comm->getSize();
    decltype( size ) root = 0;

    hmemo::HArray<PartitionId> owners;   // default is empty array

    if ( rank == root )
    {
        owners = hmemo::HArray<PartitionId>( { 0, 0, size, 0 } );   // out-of-range owner
    }

    BOOST_CHECK_THROW(
    {
        dist = generalDistributionBySingleOwners( owners, root, comm );
    }, common::Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( buildTest )
{
    using namespace hmemo;

    auto comm = Communicator::getCommunicatorPtr();
 
    auto size = comm->getSize();

    DistributionPtr dist;

    if ( size < 2 )
    {
        return;
    }

    const IndexType N = 8;

    if ( rank == 0 )
    {
        dist = generalDistribution( N, HArray<IndexType>( { 1, 3, 7, 4 } ), comm );
    }
    else if ( rank == 1 )
    {
        dist = generalDistribution( N, HArray<IndexType>( { 6, 5, 0, 2 } ), comm );
    }
    else 
    {
        dist = generalDistribution( N, HArray<IndexType>( {} ), comm );
    }
    
    // now prove that each processor has the right local values

    hmemo::HArray<PartitionId> owners;

    PartitionId root = size - 1;  // last processor

    dist->allOwners( owners, root );

    hmemo::HArray<PartitionId> expOwners( { 1, 0, 1, 0, 0, 1, 1, 0 } );

    if ( rank == root )
    {
        BOOST_TEST( hostReadAccess( owners ) == hostReadAccess( expOwners ), boost::test_tools::per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( checkTest )
{
    using namespace hmemo;

    auto comm = Communicator::getCommunicatorPtr();
 
    auto size = comm->getSize();
    auto rank = comm->getRank();

    DistributionPtr dist;

    if ( size < 2 )
    {
        BOOST_CHECK_THROW(
        {
            dist = generalDistribution( 5, HArray<IndexType>( { 0, 2, 4, 3 } ), comm );
        }, common::Exception );

        BOOST_CHECK_THROW(
        {
            dist = generalDistribution( 5, HArray<IndexType>( { 0, 2, 4, 3, 2 } ), comm );
        }, common::Exception );

        dist = generalDistribution( 5, HArray<IndexType>( { 0, 2, 4, 3, 1 } ), comm );
    }

    if ( size == 2 )
    {
        HArray<IndexType> myIndexes;

        // generate wrong data where global index 2 has two owners
        if ( rank == 0 )
        {
            myIndexes = HArray<IndexType>( { 0, 1, 2 } );
        }
        else
        {
            myIndexes = HArray<IndexType>( { 2, 3, 4 } ); 
        }

        BOOST_CHECK_THROW(
        {
            dist = generalDistribution( 5, myIndexes, comm );
        }, common::Exception );

        // generate wrong data where global index 2 has no owners
        if ( rank == 0 )
        {
            myIndexes = HArray<IndexType>( { 0, 1 } );
        }
        else
        {
            myIndexes = HArray<IndexType>( { 3, 4 } ); 
        }

        BOOST_CHECK_THROW(
        {
            dist = generalDistribution( 5, myIndexes, comm );
        }, common::Exception );

        // generate wrong data where global index 2 has two and 3 no owners
        // but we have at least sum ( myIndexes.size() ) = globalSize

        if ( rank == 0 )
        {
            myIndexes = HArray<IndexType>( { 0, 1, 2 } );
        }
        else
        {
            myIndexes = HArray<IndexType>( { 2, 4 } ); 
        }

        BOOST_CHECK_THROW(
        {
            dist = generalDistribution( 5, myIndexes, comm );
        }, common::Exception );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( uncheckedTest )
{
    using namespace hmemo;
    using namespace utilskernel;

    const IndexType N = 10;

    auto comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    // owned indexes for block( N * size )

    HArray<IndexType> myIndexes;  // { rank * N, ..., rank * N + N - 1 } 
    HArrayUtils::setSequence<IndexType>( myIndexes, N * rank, 1, N );

    // save the pointer to the allocated host data to check for correct move 

    const IndexType* ptr = hostReadAccess( myIndexes );

    auto dist = generalDistributionUnchecked( N * size, std::move( myIndexes ), comm );

    // due to move allocated memory is reused

    BOOST_CHECK_EQUAL( ptr, hostReadAccess( dist->getMyIndexes() ) );

    HArray<IndexType> globalIndexes;  // { 0, N, 2*N, (size-1) * N }
    HArrayUtils::setSequence<IndexType>( globalIndexes, 0, N, size );

    HArray<PartitionId> expectedOwners;
    HArrayUtils::setSequence<PartitionId>( expectedOwners, 0, 1, comm->getSize() );

    auto owners = dist->owner( globalIndexes );

    BOOST_TEST( hostReadAccess( expectedOwners ) == hostReadAccess( owners ), boost::test_tools::per_element() );

    if ( size < 2 ) 
    {
        return;
    }

    if ( rank < 2 )
    {
         // swap owned indexes for processors 0 and 1

         HArrayUtils::setSequence<IndexType>( myIndexes, N * ( 1 - rank ), 1, N );
         dist = generalDistributionUnchecked( N * size, std::move( myIndexes ), comm );
    }
    else
    {
         // mandatory to get all distributions in the same state, otherwise it hangs
         dist = generalDistributionUnchecked( std::move( dist ) );
    }

    // update the expected owners

    {
        auto wExpectedOwners = hostWriteAccess( expectedOwners );
        wExpectedOwners[0] = 1;
        wExpectedOwners[1] = 0;
    }

    dist->computeOwners( owners, globalIndexes );

    BOOST_TEST( hostReadAccess( expectedOwners ) == hostReadAccess( owners ), boost::test_tools::per_element() );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( uncheckedTestMove )
{
    using namespace hmemo;
    using namespace utilskernel;

    const IndexType N = 10;

    auto comm = Communicator::getCommunicatorPtr();

    PartitionId rank = comm->getRank();
    PartitionId size = comm->getSize();

    // owned indexes for block( N * size )

    HArray<IndexType> myIndexes;  // { rank * N, ..., rank * N + N - 1 } 
    HArrayUtils::setSequence<IndexType>( myIndexes, N * rank, 1, N );

    // save the pointer to the allocated host data to check for correct move 

    const IndexType* ptr = hostReadAccess( myIndexes );

    auto dist = generalDistributionUnchecked( N * size, std::move( myIndexes ), comm );

    // due to move allocated memory is reused

    BOOST_CHECK_EQUAL( ptr, hostReadAccess( dist->getMyIndexes() ) );

    auto dist1 = generalDistributionUnchecked( std::move( dist ) );

    // due to move allocated memory is reused

    BOOST_CHECK_EQUAL( ptr, hostReadAccess( dist1->getMyIndexes() ) );

    auto dist2 = generalDistributionUnchecked( dist1 );

    BOOST_CHECK_EQUAL( ptr, hostReadAccess( dist1->getMyIndexes() ) );
    BOOST_CHECK( ptr != hostReadAccess( dist2->getMyIndexes() ) );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
