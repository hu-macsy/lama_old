/**
 * @file TimerTest.cpp
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
 * @brief Test routines for the class Timer.
 * @author Matthias Makulla
 * @date 03.02.2012
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/logger/Timer.hpp>

#include <scai/logging.hpp>

#include <scai/common/macros/throw.hpp>
#include <scai/common/Walltime.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai;
using namespace solver;

using common::Exception;

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_SUITE( TimerTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.TimerTest" )

static int TIMER_ACCURACY = 10;

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
    common::Walltime::sleep( 100 );  // 100 ms, 0.1 s
    double time1 = timer.getTime( "Timer" );
    BOOST_CHECK( 0.0 < timer.getTime( "Timer" ) );
    //Call reset, but do not stop the timer
    timer.reset( "Timer" );
    common::Walltime::sleep( 10 );  // 10 ms, 0.01 s
    double time2 = timer.getTime( "Timer" );
    timer.stop( "Timer" );
    BOOST_CHECK( time2 < time1 );
    BOOST_CHECK_CLOSE( time1, 0.1, TIMER_ACCURACY );
    BOOST_CHECK_CLOSE( time2, 0.01, TIMER_ACCURACY );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( StartAndStopTimerTest )
{
    Timer timer;
    timer.start( "TestTimer" );
    common::Walltime::sleep( 100 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResumeTimerTest )
{
    Timer timer;
    timer.start( "TestTimer2" );
    timer.start( "TestTimer" );
    common::Walltime::sleep( 100 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    common::Walltime::sleep( 100 );
    timer.start( "TestTimer" );
    common::Walltime::sleep( 100 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    common::Walltime::sleep( 100 );
    timer.start( "TestTimer" );
    common::Walltime::sleep( 100 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.3, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    timer.stop( "TestTimer2" );
    BOOST_CHECK_CLOSE( 0.5, timer.getTime( "TestTimer2" ), TIMER_ACCURACY );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( ResetTimerTest )
{
    Timer timer;
    timer.start( "TestTimer" );
    common::Walltime::sleep( 100 );  // 100 ms, 0.1 s
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    timer.stopAndReset( "TestTimer" );
    timer.start( "TestTimer" );
    common::Walltime::sleep( 100 );
    timer.stop( "TestTimer" );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    timer.stopAndReset( "TestTimer" );
}

/* --------------------------------------------------------------------- */

BOOST_AUTO_TEST_CASE( GetTimeTest )
{
    Timer timer;
    timer.start( "TestTimer" );
    timer.start( "TestTimer2" );
    common::Walltime::sleep( 100 );
    BOOST_CHECK_CLOSE( 0.1, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    common::Walltime::sleep( 100 );
    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer" ), TIMER_ACCURACY );
    BOOST_CHECK_CLOSE( 0.2, timer.getTime( "TestTimer2" ), TIMER_ACCURACY );
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
