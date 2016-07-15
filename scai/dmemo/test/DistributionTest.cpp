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

    AllDistributions()
    {
        CommunicatorPtr comm = Communicator::getCommunicatorPtr();

        const IndexType globalSize = 17;

        std::vector<std::string> values;

        Distribution::getCreateValues( values );

        for ( size_t i = 0; i < values.size(); ++i )
        {
            DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

            BOOST_CHECK_EQUAL( dist->getKind(), values[i] );

            push_back( dist );
        } 
    }
};

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( localSizeTest )
{
    AllDistributions allDist;

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
    const IndexType globalSize = 17;
    std::vector<std::string> values;
    Distribution::getCreateValues( values );
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );
        SCAI_LOG_INFO( logger, *comm << ": local2GlobalTest, dist = " << *dist )

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
    const IndexType globalSize = 17;
    std::vector<std::string> values;
    Distribution::getCreateValues( values );
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

        SCAI_LOG_INFO( logger, *comm << ": global2LocalTest, dist = " << *dist )

        for ( IndexType i = 0; i < dist->getLocalSize(); i++ )
        {
            BOOST_CHECK_EQUAL( i, dist->global2local( dist->local2global( i ) ) );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ownedIndexesTest )
{
    AllDistributions allDist;

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
    AllDistributions allDist;

    for ( size_t i = 0; i < allDist.size(); ++i )
    {
        DistributionPtr dist = allDist[i];

        IndexType nGlobal = dist->getGlobalSize();

        utilskernel::LArray<PartitionId> indexes;
        utilskernel::HArrayUtils::setOrder( indexes, nGlobal );

        utilskernel::LArray<PartitionId> owners1;
        utilskernel::LArray<PartitionId> owners2;

        dist->computeOwners( owners1, indexes );                 // call the efficient derived class method
        dist->Distribution::computeOwners( owners2, indexes );   // call the straight forw from base class

        BOOST_REQUIRE_EQUAL( nGlobal, owners1.size() );
        BOOST_REQUIRE_EQUAL( nGlobal, owners2.size() );

        hmemo::ReadAccess<IndexType> rOwners1( owners1 );
        hmemo::ReadAccess<IndexType> rOwners2( owners2 );

        for ( IndexType i = 0; i < nGlobal; ++i )
        {
            BOOST_CHECK_EQUAL( rOwners1[i], rOwners2[i] );
        }
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    const IndexType globalSize = 17;
    std::vector<std::string> values;
    Distribution::getCreateValues( values );
    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );
        SCAI_LOG_INFO( logger, *comm << ": writeAt, dist = " << *dist )
        std::ostringstream out;
        out << *dist;
        BOOST_CHECK( out.str().length() > 0 );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */


