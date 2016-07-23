/**
 * @file DistributionTest.cpp
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
 * @brief Tests that will be applied to all registered distributions of the factory.
 * @author Thomas Brandes
 * @date 30.07.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/dmemo.hpp>
#include <scai/dmemo/GenBlockDistribution.hpp>
#include <scai/dmemo/GeneralDistribution.hpp>
#include <scai/utilskernel.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DistributionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DistributionTest" )

/* --------------------------------------------------------------------- */

class AllDistributions : public std::vector<DistributionPtr> 
{
public:

    AllDistributions( const IndexType globalSize )
    {
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        std::vector<std::string> values;

        Distribution::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

            BOOST_CHECK_EQUAL( dist->getKind(), values[i] );

            push_back( dist );
        } 

        utilskernel::LArray<PartitionId> owners;

        {
            PartitionId owner = 315;
            PartitionId nPartitions = comm->getSize();

            hmemo::WriteOnlyAccess<PartitionId> wOwners( owners, globalSize );

            for ( IndexType i = 0; i < globalSize; ++i )
            {
                owner = owner * 119 % 185;
                wOwners[i] = owner % nPartitions;
            }
        }

        push_back( DistributionPtr( new GeneralDistribution( owners, comm ) ) );

        float weight = static_cast<float>( comm->getRank() + 1 );

        push_back( DistributionPtr( new GenBlockDistribution( globalSize, weight, comm ) ) );
    }

private:

    AllDistributions();

};

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( localSizeTest )
{
    AllDistributions allDist( 17 );

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

BOOST_AUTO_TEST_CASE( local2GlobalTest )
{
    AllDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        SCAI_LOG_INFO( logger, "local2GlobalTest, dist = " << *dist )

        for ( IndexType i = 0; i < dist->getGlobalSize(); i++ )
        {
            if ( dist->isLocal( i ) )
            {
                BOOST_CHECK_EQUAL( i, dist->local2global( dist->global2local( i ) ) );
            }
            else
            {
                BOOST_CHECK_EQUAL( nIndex, dist->global2local( i ) );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( global2LocalTest )
{
    AllDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        const Communicator& comm = dist->getCommunicator();

        SCAI_LOG_INFO( logger, comm << ": global2LocalTest, dist = " << *dist )

        for ( IndexType i = 0; i < dist->getLocalSize(); i++ )
        {
            BOOST_CHECK_EQUAL( i, dist->global2local( dist->local2global( i ) ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ownedIndexesTest )
{
    AllDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        utilskernel::LArray<IndexType> myIndexes1;
        utilskernel::LArray<IndexType> myIndexes2;

        dist->getOwnedIndexes( myIndexes1 );                 // call it for block distribution
        dist->Distribution::getOwnedIndexes( myIndexes2 );   // call if from base class

        IndexType nLocal = dist->getLocalSize();

        BOOST_REQUIRE_EQUAL( nLocal, myIndexes1.size() );
        BOOST_REQUIRE_EQUAL( nLocal, myIndexes2.size() );

        hmemo::ReadAccess<IndexType> rIndexes1( myIndexes1 );
        hmemo::ReadAccess<IndexType> rIndexes2( myIndexes2 );

        for ( IndexType i = 0; i < nLocal; ++i )
        {
            BOOST_CHECK_EQUAL( rIndexes1[i], rIndexes2[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( computeOwnersTest )
{
    AllDistributions allDist( 17 );

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];
 
        const PartitionId rank = dist->getCommunicator().getRank();
        const PartitionId root = dist->getCommunicator().getSize() / 2;

        IndexType nGlobal = dist->getGlobalSize();

        utilskernel::LArray<PartitionId> indexes;
        utilskernel::HArrayUtils::setOrder( indexes, nGlobal );

        utilskernel::LArray<PartitionId> owners1;
        utilskernel::LArray<PartitionId> owners2;
        utilskernel::LArray<PartitionId> owners3;

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
            BOOST_REQUIRE_EQUAL( 0, owners3.size() );
        }

        hmemo::ReadAccess<IndexType> rOwners1( owners1 );
        hmemo::ReadAccess<IndexType> rOwners2( owners2 );
        hmemo::ReadAccess<IndexType> rOwners3( owners3 );

        for ( IndexType i = 0; i < nGlobal; ++i )
        {
            BOOST_CHECK_EQUAL( rOwners1[i], rOwners2[i] );

            if ( rank == root )
            {
                BOOST_CHECK_EQUAL( rOwners1[i], rOwners3[i] );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( getBlockDistributionSizeTest )
{
    // Test of getBlockDistributionSize() can be done by computeOwners on all processor
    // getBlockDistributionSize() != nIndex iff isAscending( owners( {0, ..., globalSize-1} ) )

    IndexType globalSizes[] = { 0, 1, 2, 3, 7, 16 };

    IndexType nCases = sizeof( globalSizes ) / sizeof( IndexType );

    for ( int k = 0; k < nCases; ++k )
    {
        AllDistributions allDist( globalSizes[k] );

        for ( size_t i = 0; i < allDist.size(); ++i )
        {
            DistributionPtr dist = allDist[i];

            IndexType nGlobal = dist->getGlobalSize();

            utilskernel::LArray<PartitionId> indexes;
            utilskernel::HArrayUtils::setOrder( indexes, nGlobal );

            utilskernel::LArray<PartitionId> owners;

            dist->computeOwners( owners, indexes );

            bool ascending = true;
            bool isSorted = utilskernel::HArrayUtils::isSorted( owners, ascending );

            IndexType bs = dist->getBlockDistributionSize();

            SCAI_LOG_DEBUG( logger, *dist << ", owners sorted = " << isSorted << ", bs = " << bs )

            if ( isSorted )
            {
                BOOST_CHECK_EQUAL( bs, dist->getLocalSize() );
            }
            else
            {
                BOOST_CHECK_EQUAL( bs, nIndex );
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateTest )
{
    const IndexType N = 30;  // global size 

    AllDistributions allDist( N );

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
                wLocalValues[i] = dist->local2global( i );
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
    const IndexType globalN = 11; // global size, ToDo: this test fails very strange with N = 15, 17
    const IndexType repN = 4;

    AllDistributions allDist( globalN );

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
                IndexType val = dist->local2global( i );

                for ( IndexType k = 0; k < repN; ++k )
                {
                    wLocalValues[ repN * i + k ] = val;
                }
            } 
        }

        hmemo::HArray<IndexType> allValues( repN * globalN, nIndex );

        // Now replicate the local values 

        {
            hmemo::WriteAccess<IndexType> wAllValues( allValues, repN * globalN );
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
                    SCAI_LOG_ERROR( logger, dist->getCommunicator() << ": dist = " << *dist 
                                            << ", wrong at i = " << i << " of " << globalN 
                                            << ", k = " << k << " of repN = " << repN )
                }
            }
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( replicateRaggedTest )
{
    const IndexType globalN = 15;  // global size 
    const IndexType repN = 3;

    AllDistributions allDist( globalN );

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
                IndexType val = dist->local2global( i );

                for ( IndexType k = 0; k < repN; ++k )
                {
                    wLocalValues[ repN * i + k ] = val;
                }
            } 
        }

        hmemo::HArray<IndexType> offsets( globalN, repN );
        IndexType totalValues = utilskernel::HArrayUtils::scan( offsets );

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
    AllDistributions allDist( 17 );

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
    AllDistributions allDist1( 17 );
    AllDistributions allDist2( 17 );
    AllDistributions allDist3( 18 );

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


