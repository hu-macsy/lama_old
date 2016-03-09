/**
 * @file TFQMRT.cpp
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
 * @brief TFQMRT.cpp
 * @author: Jan Ecker
 * @date 09.03.2016
 * @since 2.0.0
 **/

#include <boost/test/unit_test.hpp>

#include <scai/solver/TFQMR.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( TFQMRTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.TFQMRTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<TFQMR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    TFQMR TFQMRSolver( "TFQMRTestSolver", slogger );
    BOOST_CHECK_EQUAL( TFQMRSolver.getId(), "TFQMRTestSolver" );

    TFQMR TFQMRSolver2( "TFQMRTestSolver2" );
    BOOST_CHECK_EQUAL( TFQMRSolver2.getId(), "TFQMRTestSolver2" );

    TFQMR TFQMRSolver3( TFQMRSolver2 );
    BOOST_CHECK_EQUAL( TFQMRSolver3.getId(), "TFQMRTestSolver2" );
    BOOST_CHECK( TFQMRSolver3.getPreconditioner() == 0 );

    TFQMR TFQMRSolver4( "TFQMRSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    TFQMRSolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    TFQMRSolver4.setStoppingCriterion( criterion );

    TFQMR TFQMRSolver5( TFQMRSolver4 );
    BOOST_CHECK_EQUAL( TFQMRSolver5.getId(), TFQMRSolver4.getId() );
    BOOST_CHECK_EQUAL( TFQMRSolver5.getPreconditioner()->getId(), TFQMRSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
