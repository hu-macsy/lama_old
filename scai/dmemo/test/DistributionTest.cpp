/**
 * @file DistributionTest.cpp
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
 * @brief Tests that will be applied to all registered distributions of the factory.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/utilskernel.hpp>

#include <scai/dmemo/test/TestDistributions.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DistributionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DistributionTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( localSizeTest )
{
    TestDistributions allDist( 16 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        const Communicator& comm = dist->getCommunicator();

        SCAI_LOG_INFO( logger, comm << ": localSizeTest, dist = " << *dist )

        // Do not use comm for reductions as NoDistribution has NoCommunicator

        IndexType sumLocalSizes = dist->getCommunicator().sum( dist->getLocalSize() );

        BOOST_CHECK_EQUAL( dist->getGlobalSize(), sumLocalSizes );
    }
}

/* --------------------------------------------------------------------- */


BOOST_AUTO_TEST_CASE( maxLocalSizeTest )
{
    TestDistributions allDist( 16 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        const Communicator& comm = dist->getCommunicator();

        SCAI_LOG_INFO( logger, comm << ": maxLocalSizeTest, dist = " << *dist )

        // most distribution know an efficient implementation
        IndexType maxLocalSizeComputed = dist->getMaxLocalSize();

        // this computes the smallest maximum
        IndexType maxLocalSizeExpected = comm.max( dist->getLocalSize() );

        BOOST_CHECK_EQUAL( maxLocalSizeExpected, maxLocalSizeComputed );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( local2GlobalTest )
{
    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        SCAI_LOG_INFO( logger, "local2GlobalTest, dist = " << *dist )

        for ( IndexType i = 0; i < dist->getGlobalSize(); i++ )
        {
            if ( dist->isLocal( i ) )
            {
                BOOST_CHECK_EQUAL( i, dist->local2Global( dist->global2Local( i ) ) );
            }
            else
            {
                BOOST_CHECK_EQUAL( invalidIndex, dist->global2Local( i ) );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( global2LocalTest )
{
    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        const Communicator& comm = dist->getCommunicator();

        SCAI_LOG_INFO( logger, comm << ": global2LocalTest, dist = " << *dist )

        for ( IndexType i = 0; i < dist->getLocalSize(); i++ )
        {
            BOOST_CHECK_EQUAL( i, dist->global2Local( dist->local2Global( i ) ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( global2LocalVTest )
{
    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        const Communicator& comm = dist->getCommunicator();

        SCAI_LOG_INFO( logger, comm << ": global2LocalTest, dist = " << *dist )

        hmemo::HArray<IndexType> globalIndexes;
        utilskernel::HArrayUtils::setOrder( globalIndexes, dist->getGlobalSize() );

        hmemo::HArray<IndexType> localIndexes;
        dist->global2LocalV( localIndexes, globalIndexes );

        auto rLocal = hostReadAccess( localIndexes );

        IndexType countLocal = 0;

        for ( IndexType i = 0; i < dist->getGlobalSize(); ++i )
        {
            if ( rLocal[i] != invalidIndex )
            {
                countLocal++;
            }
        }
 
        BOOST_CHECK_EQUAL( countLocal, dist->getLocalSize() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ownedIndexesTest )
{
    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        hmemo::HArray<IndexType> myIndexes1;
        hmemo::HArray<IndexType> myIndexes2;

        dist->getOwnedIndexes( myIndexes1 );                 // call it for block distribution
        dist->Distribution::getOwnedIndexes( myIndexes2 );   // call if from base class

        IndexType nLocal = dist->getLocalSize();

        BOOST_REQUIRE_EQUAL( nLocal, myIndexes1.size() );
        BOOST_REQUIRE_EQUAL( nLocal, myIndexes2.size() );

        BOOST_TEST( hmemo::hostReadAccess( myIndexes1 ) == hmemo::hostReadAccess( myIndexes2 ),
                    boost::test_tools::per_element() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( anyAddressingTest )
{
    TestDistributions allDist( 12 );   // allows for grid 2 x 2 x 3

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        dist->enableAnyAddressing();

        IndexType nGlobal = dist->getGlobalSize();
        IndexType nP      = dist->getCommunicator().getSize();

        // own counter array to check for good local indexes on each partition

        std::unique_ptr<IndexType[]> counts( new IndexType[nP] );

        for ( IndexType iP = 0; iP < nP; ++iP )
        {
            counts[iP] = 0;
        }

        for ( IndexType globalIndex = 0; globalIndex < nGlobal; ++globalIndex )
        {
            PartitionId owner = dist->getAnyOwner( globalIndex );
            BOOST_CHECK( owner < nP );
            PartitionId localIndex = dist->getAnyLocalIndex( globalIndex, owner );
            BOOST_CHECK_EQUAL( counts[owner], localIndex );
            counts[owner]++;
            PartitionId globalIndex1 = dist->getAnyGlobalIndex( localIndex, owner );
            BOOST_CHECK_EQUAL( globalIndex, globalIndex1 );
        }

        for ( IndexType iP = 0; iP < nP; ++iP )
        {
            BOOST_CHECK_EQUAL( counts[iP], dist->getAnyLocalSize( iP ) );
        }

        // now check for correct permutations

        hmemo::HArray<IndexType> offsets;
        hmemo::HArray<IndexType> perm;

        dist->getAnyGlobal2Local( offsets, perm );

        BOOST_REQUIRE_EQUAL( nP + 1, offsets.size() );
        BOOST_REQUIRE_EQUAL( nGlobal, perm.size() );

        for ( IndexType iP = 0; iP < nP; ++iP )
        {
            BOOST_CHECK_EQUAL( counts[iP], offsets[iP + 1] - offsets[iP] );
        }

        for ( IndexType globalIndex = 0; globalIndex < nGlobal; ++globalIndex )
        {
            PartitionId owner = dist->getAnyOwner( globalIndex );
            IndexType localIndex1 = dist->getAnyLocalIndex( globalIndex, owner );
            IndexType localIndex2 = perm[ globalIndex ] - offsets[owner];
            SCAI_LOG_TRACE( logger, "Global index = " << globalIndex << ", owner = " << owner
                            << ", local1 = " << localIndex1 << ", local2 = " << localIndex2 )
            BOOST_CHECK_EQUAL( localIndex1, localIndex2 );
        }

        BOOST_REQUIRE_EQUAL( nP + 1, offsets.size() );
        BOOST_REQUIRE_EQUAL( nGlobal, perm.size() );

        dist->getAnyLocal2Global( offsets, perm );

        for ( IndexType iP = 0; iP < nP; ++iP )
        {
            BOOST_CHECK_EQUAL( counts[iP], offsets[iP + 1] - offsets[iP] );
        }

        for ( IndexType iP = 0; iP < nP; ++iP )
        {
            IndexType nLocal = dist->getAnyLocalSize( iP );
            IndexType offset = offsets[iP];

            for ( IndexType localIndex = 0; localIndex < nLocal; ++localIndex )
            {
                IndexType globalIndex1 = dist->getAnyGlobalIndex( localIndex, iP );
                IndexType globalIndex2 = perm[localIndex + offset];
                SCAI_LOG_TRACE( logger, "iP = " << iP << ", local index = " << localIndex
                                << ", global1 = " << globalIndex1 << ", global2 = " << globalIndex2 )
                BOOST_CHECK_EQUAL( globalIndex1, globalIndex2 );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( findOwnerTest )
{
    const IndexType globalK = 6;   // element for which we find owner

    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        const auto& dist = *allDist[i];

        // take the expected owner from the array version

        auto owners = dist.owner( hmemo::HArray<IndexType>( { globalK } ) );
        PartitionId expectedOwner = owners[0];

        PartitionId owner = dist.findOwner( globalK );
        BOOST_CHECK_EQUAL( expectedOwner, owner );

        // force the use of the default method of base class
        owner = dist.Distribution::findOwner( globalK );
        BOOST_CHECK_EQUAL( expectedOwner, owner );

        if ( dist.hasAnyAddressing() )
        {
            owner = dist.getAnyOwner( globalK );
            BOOST_CHECK_EQUAL( expectedOwner, owner );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( computeOwnersTest )
{
    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        const PartitionId rank = dist->getCommunicator().getRank();
        const PartitionId root = dist->getCommunicator().getSize() / 2;

        IndexType nGlobal = dist->getGlobalSize();

        hmemo::HArray<IndexType> indexes;
        utilskernel::HArrayUtils::setOrder( indexes, nGlobal );

        hmemo::HArray<PartitionId> owners1;
        hmemo::HArray<PartitionId> owners2;
        hmemo::HArray<PartitionId> owners3;

        dist->computeOwners( owners1, indexes );                 // call the efficient derived class method
        dist->Distribution::computeOwners( owners2, indexes );   // call the straight forw from base class
        dist->allOwners( owners3, root );                        // allOwners, result only @ root

        BOOST_REQUIRE_EQUAL( nGlobal, owners1.size() );
        BOOST_REQUIRE_EQUAL( nGlobal, owners2.size() );

        if ( rank == root )
        {
            // only root processor gets the owners via allOwners

            BOOST_REQUIRE_EQUAL( nGlobal, owners3.size() );
        }
        else
        {
            BOOST_REQUIRE_EQUAL( IndexType( 0 ), owners3.size() );
        }

        BOOST_TEST( hostReadAccess( owners1 ) == hostReadAccess( owners2 ), boost::test_tools::per_element() );

        if ( rank == root )
        {
            BOOST_TEST( hostReadAccess( owners1 ) == hostReadAccess( owners3 ), boost::test_tools::per_element() );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getBlockDistributionSizeTest )
{
    // Test of getBlockDistributionSize() can be done by computeOwners on all processor
    // getBlockDistributionSize() != invalidIndex iff isAscending( owners( {0, ..., globalSize-1} ) )

    IndexType globalSizes[] = { 0, 1, 2, 3, 7, 16 };

    int nCases = sizeof( globalSizes ) / sizeof( IndexType );

    for ( int k = 0; k < nCases; ++k )
    {
        TestDistributions allDist( globalSizes[k] );

        for ( size_t i = 0; i < allDist.size(); ++i )
        {
            DistributionPtr dist = allDist[i];

            IndexType nGlobal = dist->getGlobalSize();

            hmemo::HArray<IndexType> indexes;
            utilskernel::HArrayUtils::setOrder( indexes, nGlobal );

            hmemo::HArray<PartitionId> owners;

            dist->computeOwners( owners, indexes );

            // check for sorted owners, e.g. 0 0 0 1 1 1 1 2 2 2 2 3 3 3 indicates block dist

            bool isSorted = utilskernel::HArrayUtils::isSorted( owners, common::CompareOp::LE );

            IndexType bs = dist->getBlockDistributionSize();

            SCAI_LOG_DEBUG( logger, *dist << ", owners sorted = " << isSorted << ", bs = " << bs )

            if ( isSorted )
            {
                if ( bs == invalidIndex )
                {
                    // might happen for grid distributions
                    SCAI_LOG_WARN( logger, "Owners sorted, but not block distribution: " << *dist )
                }
                else
                {
                    BOOST_CHECK_EQUAL( bs, dist->getLocalSize() );
                }
            }
            else
            {
                BOOST_CHECK_EQUAL( bs, invalidIndex );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( convertTest )
{
    const IndexType N = 30;  // global size

    TestDistributions allDist( N );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        DistributionPtr blockedDist = dist->toBlockDistribution( comm );

        BOOST_CHECK_EQUAL( blockedDist->getGlobalSize(), dist->getGlobalSize() );
        BOOST_CHECK( blockedDist->isBlockDistributed( comm ) );

        DistributionPtr singleDist = dist->toSingleDistribution( comm );

        BOOST_CHECK_EQUAL( singleDist->getGlobalSize(), dist->getGlobalSize() );

        if ( !singleDist->isSingleDistributed( comm ) )
        {
            SCAI_LOG_ERROR( logger, "dist = " << *dist << ", single = " << *singleDist )
        }

        BOOST_CHECK( singleDist->isSingleDistributed( comm ) );

        DistributionPtr repDist = dist->toReplicatedDistribution();

        BOOST_CHECK_EQUAL( repDist->getGlobalSize(), dist->getGlobalSize() );

        if ( !repDist->isReplicated() )
        {
            SCAI_LOG_ERROR( logger, "dist = " << *dist << ", rep = " << *repDist )
        }

        BOOST_CHECK( repDist->isReplicated() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateTest )
{
    const IndexType N = 30;  // global size

    TestDistributions allDist( N );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        IndexType localN = dist->getLocalSize();

        hmemo::HArray<IndexType> allValues;
        hmemo::HArray<IndexType> localValues;

        // for test each processor fills its local values with its owned global indexes
        {
            hmemo::WriteOnlyAccess<IndexType> wLocalValues( localValues, localN );

            for ( IndexType i = 0; i < localN; ++i )
            {
                wLocalValues[i] = dist->local2Global( i );
            }
        }

        // Now replicate the local values

        {
            hmemo::WriteOnlyAccess<IndexType> wAllValues( allValues, N );
            hmemo::ReadAccess<IndexType> rLocalValues( localValues );
            dist->replicate( wAllValues.get(), rLocalValues.get() );
        }

        // the replicated array must now contain all global indexes in correct order

        hmemo::ReadAccess<IndexType> rAllValues( allValues );

        for ( IndexType i = 0; i < N; ++i )
        {
            BOOST_CHECK_EQUAL( i, rAllValues[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateNTest )
{
    const IndexType globalN = 17; // global size
    const IndexType repN = 4;

    TestDistributions allDist( globalN );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        IndexType localN = dist->getLocalSize();

        hmemo::HArray<IndexType> localValues;

        // for test each processor fills its local values with its owned global indexes
        {
            hmemo::WriteOnlyAccess<IndexType> wLocalValues( localValues, localN * repN );

            for ( IndexType i = 0; i < localN; ++i )
            {
                IndexType val = dist->local2Global( i );

                for ( IndexType k = 0; k < repN; ++k )
                {
                    wLocalValues[ repN * i + k ] = val;
                }
            }
        }

        hmemo::HArray<IndexType> allValues( repN * globalN, invalidIndex );

        // Now replicate the local values

        {
            hmemo::WriteAccess<IndexType> wAllValues( allValues );
            hmemo::ReadAccess<IndexType> rLocalValues( localValues );
            dist->replicateN( wAllValues.get(), rLocalValues.get(), repN );
        }

        // the replicated array must now contain all global indexes in correct order

        hmemo::ReadAccess<IndexType> rAllValues( allValues );

        for ( IndexType i = 0; i < globalN; ++i )
        {
            for ( IndexType k = 0; k < repN; ++k )
            {
                BOOST_CHECK_EQUAL( i, rAllValues[ repN * i + k ] );

                if ( i != rAllValues[ repN * i + k ] )
                {
                    SCAI_LOG_ERROR( logger, dist->getCommunicator() << ": dist = " << *dist <<
                                    ", wrong at i = " << i << " of " << globalN <<
                                    ", k = " << k << " of repN = " << repN <<
                                    ", rAllValues [ " << repN * i + k << " ] = " << rAllValues[ repN * i + k ] )

                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateRaggedTest )
{
    const IndexType globalN = 17;  // global size
    const IndexType repN = 4;

    TestDistributions allDist( globalN );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        IndexType localN = dist->getLocalSize();

        hmemo::HArray<IndexType> localValues;

        // for test each processor fills its local values with its owned global indexes
        {
            hmemo::WriteOnlyAccess<IndexType> wLocalValues( localValues, localN * repN );

            for ( IndexType i = 0; i < localN; ++i )
            {
                IndexType val = dist->local2Global( i );

                for ( IndexType k = 0; k < repN; ++k )
                {
                    wLocalValues[ repN * i + k ] = val;
                }
            }
        }

        hmemo::HArray<IndexType> offsets( globalN, repN );
        IndexType totalValues = utilskernel::HArrayUtils::scan1( offsets );

        hmemo::HArray<IndexType> allValues;  // result array for replicateRagged

        // Now replicate the local values

        {
            hmemo::WriteOnlyAccess<IndexType> wAllValues( allValues, repN * globalN );
            hmemo::ReadAccess<IndexType> rLocalValues( localValues );
            hmemo::ReadAccess<IndexType> rOffsets( offsets );
            dist->replicateRagged( wAllValues.get(), rLocalValues.get(), rOffsets.get() );
        }

        BOOST_CHECK_EQUAL( allValues.size(), totalValues );

        // the replicated array must now contain all global indexes in correct order

        hmemo::ReadAccess<IndexType> rAllValues( allValues );

        for ( IndexType i = 0; i < globalN; ++i )
        {
            for ( IndexType k = 0; k < repN; ++k )
            {
                BOOST_CHECK_EQUAL( i, rAllValues[ repN * i + k ] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    TestDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        SCAI_LOG_INFO( logger, "writeAt, dist = " << *dist )
        std::ostringstream out;
        out << *dist;
        BOOST_CHECK( out.str().length() > 0 );

        // verify that a derived distribution class has overridden the
        // default implementation of the base class Distriution

        std::ostringstream outBaseClass;
        dist->Distribution::writeAt( outBaseClass );
        BOOST_CHECK( out.str() != outBaseClass.str() );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( equalTest )
{
    TestDistributions allDist1( 17 );
    TestDistributions allDist2( 17 );
    TestDistributions allDist3( 18 );

    for ( size_t i = 0; i < allDist1.size(); ++i )
    {
        DistributionPtr dist1 = allDist1[i];
        DistributionPtr dist2 = allDist2[i];
        DistributionPtr dist3 = allDist3[i];

        BOOST_CHECK_EQUAL( *dist1, *dist1 );  // pointer equality
        BOOST_CHECK_EQUAL( *dist1, *dist2 );  // same distibution
        BOOST_CHECK( *dist1 != *dist3 );      // must be different due to other global size
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */


