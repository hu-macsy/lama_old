/**
 * @file CGNRTest.cpp
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
 * @brief CGNRTest.cpp
 * @author Jan Ecker
 * @date 09.03.2016
 * @since 2.0.0
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/CGNR.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( CGNRTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGNRTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<CGNR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGNR CGNRSolver( "CGNRTestSolver", slogger );
    BOOST_CHECK_EQUAL( CGNRSolver.getId(), "CGNRTestSolver" );

    CGNR CGNRSolver2( "CGNRTestSolver2" );
    BOOST_CHECK_EQUAL( CGNRSolver2.getId(), "CGNRTestSolver2" );

    CGNR CGNRSolver3( CGNRSolver2 );
    BOOST_CHECK_EQUAL( CGNRSolver3.getId(), "CGNRTestSolver2" );
    BOOST_CHECK( CGNRSolver3.getPreconditioner() == 0 );

    CGNR CGNRSolver4( "CGNRSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGNRSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    CGNRSolver4.setStoppingCriterion( criterion );

    CGNR CGNRSolver5( CGNRSolver4 );
    BOOST_CHECK_EQUAL( CGNRSolver5.getId(), CGNRSolver4.getId() );
    BOOST_CHECK_EQUAL( CGNRSolver5.getPreconditioner()->getId(), CGNRSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
