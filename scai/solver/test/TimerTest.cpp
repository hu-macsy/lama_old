/**
 * @file TimerTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief Contains the implementation of the class TimerTest.
 * @author Alexander Büchel, Matthias Makulla
 * @date 03.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/logger/Timer.hpp>
 
#include <scai/common/macros/throw.hpp>

#include <scai/solver/test/TestMacros.hpp>

//Adding support for Timers under Windows
#ifdef _WIN32
#include <windows.h>
inline void usleep( int t )
{
    Sleep( t / 1000 );
}
#elif _WIN32
#include <unistd.h>
#endif

using namespace scai::solver;
using namespace scai::hmemo;
using scai::common::Exception;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TimerTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.TimerTest" )

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    Timer timer;
    timer.start( "Solver" );
    // cannot start timer twice
    SCAI_CHECK_THROW( { timer.start( "Solver" ); }, Exception );
    // cannot stop non-existing timer
    SCAI_CHECK_THROW( { timer.stop( "Solver1" ); }, Exception );
    timer.stop( "Solver" );
    // cannot stop timer twice
    SCAI_CHECK_THROW( { timer.stop( "Solver" ); }, Exception );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResetTest )
{
    Timer timer;
    SCAI_CHECK_THROW( { timer.reset( "Timer" ) ; }, Exception );
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
    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer" ), 2 );
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
    BOOST_CHECK_THROW( timer.stopAndReset( "TestTimer" ), Exception );
    Timer timer1;
    BOOST_CHECK_THROW( timer1.reset( "TestTimer" ), Exception );
}
/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE_END();
