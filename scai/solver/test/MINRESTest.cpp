/**
 * @file MINRESTest.cpp
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
 * @brief MINRESTest.cpp
 * @author: Jan Ecker
 * @date 09.03.2016
 * @since 2.0.0
 **/

#include <boost/test/unit_test.hpp>

#include <scai/solver/MINRES.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( MINRESTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.MINRESTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<MINRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    MINRES MINRESSolver( "MINRESSolver", slogger );
    BOOST_CHECK_EQUAL( MINRESSolver.getId(), "MINRESSolver" );

    MINRES MINRESSolver2( "MINRESSolver2" );
    BOOST_CHECK_EQUAL( MINRESSolver2.getId(), "MINRESSolver2" );

    MINRES MINRESSolver3( MINRESSolver2 );
    BOOST_CHECK_EQUAL( MINRESSolver3.getId(), "MINRESSolver2" );
    BOOST_CHECK( MINRESSolver3.getPreconditioner() == 0 );

    MINRES MINRESSolver4( "MINRESSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    MINRESSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    MINRESSolver4.setStoppingCriterion( criterion );

    MINRES MINRESSolver5( MINRESSolver4 );
    BOOST_CHECK_EQUAL( MINRESSolver5.getId(), MINRESSolver4.getId() );
    BOOST_CHECK_EQUAL( MINRESSolver5.getPreconditioner()->getId(), MINRESSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
