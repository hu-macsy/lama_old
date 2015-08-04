/**
 * @file TaskSyncTokenTest.hpp
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
 * @brief Contains the implementation of the class TaskSyncTokenTest.
 * @author Alexander Büchel, Thomas Brandes
 * @date 02.02.2012
 * @since 1.0.0
 */

#include <boost/test/unit_test.hpp>

#include <tasking/TaskSyncToken.hpp>

#include <test/TestMacros.hpp>
#include <common/bind.hpp>

using namespace tasking;
using namespace common;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TaskSyncTokenTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.TaskSyncTokenTest" )

/* ----------------------------------------------------------------------- */

void threadMethod( const int in, int& out )
{
    out = in;
}

/* ----------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( runTest )
{
    int in = 1;
    int out = 0;
    //TODO: Check /*this*/
    TaskSyncToken testToken( bind( &TaskSyncTokenTest::threadMethod, /*this,*/in, ref( out ) ) );
    testToken.wait();
    BOOST_CHECK_EQUAL( in, out );
}

/* ----------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( writeAtTest )
{
    int in = 1;
    int out = 0;
    TaskSyncToken testToken( bind( &TaskSyncTokenTest::threadMethod, /*this,*/in, ref( out ) ) );
    LAMA_WRITEAT_TEST( testToken );
}
/* ----------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
