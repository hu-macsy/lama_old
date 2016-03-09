/**
 * @file RichardsonTest.cpp
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
 * @brief Contains the implementation of the class RichardsonTest.
 * @author: 
 * @date 17.04.2015
 * @since 
 **/
#include <boost/test/unit_test.hpp>

#include <scai/solver/Richardson.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( RichardsonTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.RichardsonTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<Richardson>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    Richardson rSolver( "RichardsonTestSolver", slogger );
    BOOST_CHECK_EQUAL( rSolver.getId(), "RichardsonTestSolver" );

    Richardson rSolver2( "RichardsonTestSolver2" );
    BOOST_CHECK_EQUAL( rSolver2.getId(), "RichardsonTestSolver2" );

    Richardson rSolver3( rSolver2 );
    BOOST_CHECK_EQUAL( rSolver3.getId(), "RichardsonTestSolver2" );
    BOOST_CHECK( rSolver3.getPreconditioner() == 0 );

    Richardson rSolver4( "RichardsonTestSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    rSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    rSolver4.setStoppingCriterion( criterion );

    Richardson rSolver5( rSolver4 );
    BOOST_CHECK_EQUAL( rSolver5.getId(), rSolver4.getId() );
    BOOST_CHECK_EQUAL( rSolver5.getPreconditioner()->getId(), rSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
