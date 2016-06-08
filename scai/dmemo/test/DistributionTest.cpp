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

using namespace scai::dmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( DistributionTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.DistributionTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( localSizeTest )
{
    const IndexType globalSize = 17;

    std::vector<std::string> values;

    Distribution::getCreateValues( values );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

        BOOST_CHECK_EQUAL( dist->getKind(), values[i] );

        SCAI_LOG_INFO( logger, *comm << ": localSizeTest, dist = " << *dist )

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

BOOST_AUTO_TEST_CASE( printDistributionVector )
{
    const IndexType globalSize = 17;

    std::vector<std::string> values;

    Distribution::getCreateValues( values );

    CommunicatorPtr comm = Communicator::getCommunicatorPtr();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistributionPtr( values[i], comm, globalSize ) );

        // ToDo: does not test the content of these files

        std::string fileName = "distribution.";
        fileName += dist->getKind();

        SCAI_LOG_INFO( logger, *comm << ": printDistributionVector, dist = " << *dist << ", filename = " << fileName )

        dist->printDistributionVector( fileName );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */


