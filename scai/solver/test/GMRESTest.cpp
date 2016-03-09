/**
 * @file GMRESTest.cpp
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
 * @brief Contains the implementation of the class GMRESTest.
 * @author: Malte Förster
 * @date 10.04.2012
 * @since 1.0.0
 **/

#include <boost/test/unit_test.hpp>

#include <scai/solver/test/TestMacros.hpp>

#include <scai/solver/GMRES.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( GMRESTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.GMRESTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<GMRES>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    GMRES GMRESSolver( "GMRESSolver", slogger );
    BOOST_CHECK_EQUAL( GMRESSolver.getId(), "GMRESSolver" );

    GMRES GMRESSolver2( "GMRESSolver2" );
    BOOST_CHECK_EQUAL( GMRESSolver2.getId(), "GMRESSolver2" );

    GMRES GMRESSolver3( GMRESSolver2 );
    BOOST_CHECK_EQUAL( GMRESSolver3.getId(), "GMRESSolver2" );
    BOOST_CHECK( GMRESSolver3.getPreconditioner() == 0 );

    GMRES GMRESSolver4( "GMRESSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    GMRESSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    GMRESSolver4.setStoppingCriterion( criterion );

    GMRES GMRESSolver5( GMRESSolver4 );
    BOOST_CHECK_EQUAL( GMRESSolver5.getId(), GMRESSolver4.getId() );
    BOOST_CHECK_EQUAL( GMRESSolver5.getPreconditioner()->getId(), GMRESSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
