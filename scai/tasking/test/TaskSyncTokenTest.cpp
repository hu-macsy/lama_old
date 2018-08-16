/**
 * @file TaskSyncTokenTest.cpp
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
 * @brief Extensive test program for TaskSyncToken
 * @author Thomas Brandes
 * @date 11.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <scai/common/test/TestMacros.hpp>

#include <memory>
#include <functional>

using std::bind;

using namespace scai::common;
using namespace scai::tasking;

using std::bind;
using std::ref;

static const int WORKLOAD = 1000000;

/* ----------------------------------------------------------------------- */

static void work( int& out, const int in )
{
    out = in - 1;
    int factor = in % 4 + 1;

    // just do some stupid work, workload depends on in

    for ( int i = 0; i < WORKLOAD * factor; i++ )
    {
        int dir = i & 1;

        if ( dir )
        {
            out += 13;
        }
        else
        {
            out -= 13;
        }
    }

    out = out + 1;
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TaskSyncTokenTest )

/* --------------------------------------------------------------------- */

SCAI_LOG_DEF_LOGGER( logger, "Test.TaskSyncTokenTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( constructorTest )
{
    int in  = 15;
    int out = 3;
    {
        TaskSyncToken token( bind( &work, std::ref( out ), in ) );
    }
    BOOST_CHECK_EQUAL( in, out );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( runTest )
{
    for ( int i = 0; i < 10; ++i )
    {
        int out = 3;
        TaskSyncToken token( bind( &work, std::ref( out ), i ) );
        token.wait();
        BOOST_CHECK_EQUAL( i, out );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    int out = 0;
    TaskSyncToken testToken( bind( &work, std::ref( out ), 1 ) );
    SCAI_COMMON_WRITEAT_TEST( testToken );
    testToken.wait();
    BOOST_CHECK_EQUAL( 1, out );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( fullTest )
{
    const int N = 15;
    SyncToken* tokenArray[15];
    int out[15];

    for ( int i = 0; i < N; ++i )
    {
        tokenArray[i] = new TaskSyncToken( std::bind( &work, std::ref( out[i] ), i ) );
    }

    for ( int i = 0; i < N; ++i )
    {
        delete tokenArray[i];
        BOOST_CHECK_EQUAL( i, out[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
