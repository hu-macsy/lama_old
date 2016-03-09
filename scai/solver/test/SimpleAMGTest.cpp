/**
 * @file SimpleAMGTest.cpp
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
 * @brief Contains the implementation of the class SimpleAMGTest.cpp
 * @author: Alexander Büchel, Robin Rehrmann
 * @date 22.02.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>

#include <scai/solver/SimpleAMG.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( SimpleAMGTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.SimpleAMGTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<SimpleAMG>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    SimpleAMG SimpleAMGSolver( "SimpleAMGSolver", slogger );
    BOOST_CHECK_EQUAL( SimpleAMGSolver.getId(), "SimpleAMGSolver" );

    SimpleAMG SimpleAMGSolver2( "SimpleAMGSolver2" );
    BOOST_CHECK_EQUAL( SimpleAMGSolver2.getId(), "SimpleAMGSolver2" );

    SimpleAMG SimpleAMGSolver3( SimpleAMGSolver2 );
    BOOST_CHECK_EQUAL( SimpleAMGSolver3.getId(), "SimpleAMGSolver2" );
    BOOST_CHECK( SimpleAMGSolver3.getPreconditioner() == 0 );

    SimpleAMG SimpleAMGSolver4( "SimpleAMGSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    SimpleAMGSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    SimpleAMGSolver4.setStoppingCriterion( criterion );

    SimpleAMG SimpleAMGSolver5( SimpleAMGSolver4 );
    BOOST_CHECK_EQUAL( SimpleAMGSolver5.getId(), SimpleAMGSolver4.getId() );
    BOOST_CHECK_EQUAL( SimpleAMGSolver5.getPreconditioner()->getId(), SimpleAMGSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
