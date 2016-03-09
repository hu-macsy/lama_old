/**
 * @file CGTest.cpp
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
 * @brief Contains the implementation of the class CGTest.
 * @author: Alexander Büchel, Thomas Brandes
 * @date 21.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>

#include <scai/solver/CG.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( CGTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<CG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CG cgSolver( "CGTestSolver", slogger );
    BOOST_CHECK_EQUAL( cgSolver.getId(), "CGTestSolver" );

    CG cgSolver2( "CGTestSolver2" );
    BOOST_CHECK_EQUAL( cgSolver2.getId(), "CGTestSolver2" );

    CG cgSolver3( cgSolver2 );
    BOOST_CHECK_EQUAL( cgSolver3.getId(), "CGTestSolver2" );
    BOOST_CHECK( cgSolver3.getPreconditioner() == 0 );

    CG cgSolver4( "cgSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    cgSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    cgSolver4.setStoppingCriterion( criterion );

    CG cgSolver5( cgSolver4 );
    BOOST_CHECK_EQUAL( cgSolver5.getId(), cgSolver4.getId() );
    BOOST_CHECK_EQUAL( cgSolver5.getPreconditioner()->getId(), cgSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
