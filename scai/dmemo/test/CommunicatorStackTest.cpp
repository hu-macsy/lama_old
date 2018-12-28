/**
 * @file CommunicatorStackTest.cpp
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
 * @brief Test correct working of communicator stack
 * @author Thomas Brandes
 * @date 21.12.2018
 */

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <scai/dmemo/CommunicatorStack.hpp>

using namespace scai;
using namespace dmemo;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( CommunicatorStackTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.CommunicatorTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( pushTest )
{
    CommunicatorPtr allComm = Communicator::getCommunicatorPtr();

    BOOST_REQUIRE( allComm );
    SCAI_LOG_INFO( logger, "all world comm = " << *allComm )

    PartitionId color = allComm->getRank() % 2 == 0 ? 0 : 1;
    PartitionId key   = allComm->getRank();  
 
    CommunicatorPtr comm = allComm->split( color, key );

    {
        SCAI_DMEMO_TASK( comm );
        CommunicatorPtr currComm = Communicator::getCommunicatorPtr();
        BOOST_CHECK_EQUAL( *comm, *currComm );
    }

    CommunicatorPtr currComm = Communicator::getCommunicatorPtr();
    BOOST_CHECK_EQUAL( *allComm, *currComm );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();

/* --------------------------------------------------------------------- */
