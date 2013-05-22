/**
 * @file TimerTest.cpp
 *
 * @license
 * Copyright (c) 2009-2013
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
 * @brief Contains the implementation of the class TimerTest.
 * @author: Alexander Büchel, Matthias Makulla
 * @date 03.02.2012
 * $Id$
 **/

#include <boost/test/unit_test.hpp>

#include <lama/solver/logger/Timer.hpp>
#include <lama/exception/Exception.hpp>

#include "TestMacros.hpp"

//Adding support for Timers under Windows
#ifdef _WIN32
#include <windows.h>
inline void usleep(int t)
{
    Sleep(t/1000);
}
#elif _WIN32
#include <unistd.h>
#endif

using namespace boost;
using namespace lama;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TimerTest )

LAMA_LOG_DEF_LOGGER( logger, "Test.TimerTest" );

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    Timer timer;
    timer.start( "Solver" );

    // cannot start timer twice

    LAMA_CHECK_THROW( { timer.start( "Solver" ); }, Exception );

    // cannot stop non-existing timer
    LAMA_CHECK_THROW( { timer.stop( "Solver1" ); }, Exception );

    timer.stop( "Solver" );

    // cannot stop timer twice
    LAMA_CHECK_THROW( { timer.stop( "Solver" ); }, Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResetTest )
{
    Timer timer;
    LAMA_CHECK_THROW( { timer.reset("Timer") ; }, Exception );

    timer.start( "Timer" );
    usleep( 100000 );

    double time1 = timer.getTime( "Timer" );
    BOOST_CHECK( 0.0 < timer.getTime( "Timer" ) );

    //Call reset, but do not stop the timer
    timer.reset( "Timer" );
    usleep( 10000 );

    double time2 = timer.getTime( "Timer" );
    timer.stop( "Timer" );

    BOOST_CHECK( time2 < time1 );
    BOOST_CHECK_CLOSE( time1, 0.1, 2 );
    BOOST_CHECK_CLOSE( time2, 0.01, 2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( StartAndStopTimerTest )
{
    Timer timer;

    timer.start( "TestTimer" );
    usleep( 100000 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), 2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResumeTimerTest )
{
    Timer timer;

    timer.start( "TestTimer2" );

    timer.start( "TestTimer" );
    usleep( 100000 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), 2 );

    usleep( 100000 );
    timer.start( "TestTimer" );
    usleep( 100000 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer" ), 2 );

    usleep( 100000 );
    timer.start( "TestTimer" );
    usleep( 100000 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.3, timer.getTime( "TestTimer" ), 2 );

    timer.stop( "TestTimer2" );
    BOOST_CHECK_CLOSE( 0.5, timer.getTime( "TestTimer2" ), 2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResetTimerTest )
{
    Timer timer;

    timer.start( "TestTimer" );
    usleep( 100000 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), 2 );
    timer.stopAndReset( "TestTimer" );

    timer.start( "TestTimer" );
    usleep( 100000 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), 2 );
    timer.stopAndReset( "TestTimer" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetTimeTest )
{
    Timer timer;

    timer.start( "TestTimer" );
    timer.start( "TestTimer2" );

    usleep( 100000 );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), 2 );

    usleep( 100000 );

    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer"), 2 );

    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer2" ), 2 );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( StartExceptionsTest )
{
    Timer timer;
    timer.start( "TestTimer" );
    BOOST_CHECK_THROW( timer.start( "TestTimer" ), Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( StopExceptionTest )
{
    Timer timer;
    BOOST_CHECK_THROW( timer.stop( "TestTimer" ), Exception );
    timer.start( "TestTimer" );
    timer.stop( "TestTimer" );
    BOOST_CHECK_THROW( timer.stop( "TestTimer" ), Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetTimeExceptionTest )
{
    Timer timer;
    BOOST_CHECK_THROW( timer.getTime( "TestTimer" ), Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResetExceptionTest )
{
    Timer timer;
    BOOST_CHECK_THROW( timer.stopAndReset("TestTimer"), Exception );

    Timer timer1;
    BOOST_CHECK_THROW( timer1.reset("TestTimer"), Exception );
}
/* --------------------------------------------------------------------- */BOOST_AUTO_TEST_SUITE_END();
