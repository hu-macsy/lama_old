/**
 * @file DistributionTest.cpp
 *
 * @license
 * Copyright (c) 2009-2015
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 * @endlicense
 *
 * @brief Tests that will be applied to all registered distributions of the factory.
 * @author Thomas Brandes
 * @date 30.07.2012
 * @since 1.0.0
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

    CommunicatorPtr comm = Communicator::getCommunicator();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistribution( values[i], comm, globalSize ) );

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

    CommunicatorPtr comm = Communicator::getCommunicator();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistribution( values[i], comm, globalSize ) );

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

    CommunicatorPtr comm = Communicator::getCommunicator();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistribution( values[i], comm, globalSize ) );

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

    CommunicatorPtr comm = Communicator::getCommunicator();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistribution( values[i], comm, globalSize ) );

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

    CommunicatorPtr comm = Communicator::getCommunicator();

    for ( size_t i = 0; i < values.size(); ++i )
    {
        DistributionPtr dist( Distribution::getDistribution( values[i], comm, globalSize ) );

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


