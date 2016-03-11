/**
 * @file TaskSyncTokenTest.cpp
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
 * @brief Extensive test program for TaskSyncToken
 * @author Thomas Brandes
 * @date 11.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/tasking/TaskSyncToken.hpp>

#include <memory>
#include <scai/common/bind.hpp>

using namespace scai::common;
using namespace scai::tasking;

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
        TaskSyncToken token( bind( &work, ref( out ), in ) );
    }

    BOOST_CHECK_EQUAL( in, out );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( runTest )
{
    for ( int i = 0; i < 10; ++i )
    {
        int out = 3;

        TaskSyncToken token( bind( &work, ref( out ), i ) );

        token.wait();

        BOOST_CHECK_EQUAL( i, out );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( fullTest )
{
    const int N = 15;

    SyncToken* tokenArray[15];

    int out[15];

    for ( int i = 0; i < N; ++i )
    {
        tokenArray[i] = new TaskSyncToken( bind( &work, ref( out[i] ), i ) );
    }
 
    for ( int i = 0; i < N; ++i )
    {
        delete tokenArray[i];
        BOOST_CHECK_EQUAL( i, out[i] );
    }
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
