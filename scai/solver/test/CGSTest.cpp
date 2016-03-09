/**
 * @file CGSTest.cpp
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
 * @brief CGSTest.cpp
 * @author Jan Ecker
 * @date 09.03.2016
 * @since 2.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/CGS.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( CGSTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGSTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<CGS>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGS CGSSolver( "CGSTestSolver", slogger );
    BOOST_CHECK_EQUAL( CGSSolver.getId(), "CGSTestSolver" );

    CGS CGSSolver2( "CGSTestSolver2" );
    BOOST_CHECK_EQUAL( CGSSolver2.getId(), "CGSTestSolver2" );

    CGS CGSSolver3( CGSSolver2 );
    BOOST_CHECK_EQUAL( CGSSolver3.getId(), "CGSTestSolver2" );
    BOOST_CHECK( CGSSolver3.getPreconditioner() == 0 );

    CGS CGSSolver4( "CGSSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGSSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    CGSSolver4.setStoppingCriterion( criterion );

    CGS CGSSolver5( CGSSolver4 );
    BOOST_CHECK_EQUAL( CGSSolver5.getId(), CGSSolver4.getId() );
    BOOST_CHECK_EQUAL( CGSSolver5.getPreconditioner()->getId(), CGSSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
