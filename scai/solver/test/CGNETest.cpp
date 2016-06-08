/**
 * @file CGNETest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the SCAI framework LAMA.
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
 *
 * Other Usage
 * Alternatively, this file may be used in accordance with the terms and
 * conditions contained in a signed written agreement between you and
 * Fraunhofer SCAI. Please contact our distributor via info[at]scapos.com.
 * @endlicense
 *
 * @brief CGNETest.cpp
 * @author Jan Ecker
 * @date 09.03.2016
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/CGNE.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/criteria/ResidualThreshold.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( CGNETest )

SCAI_LOG_DEF_LOGGER( logger, "Test.CGNETest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<CGNE>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );

    CGNE CGNESolver( "CGNETestSolver", slogger );
    BOOST_CHECK_EQUAL( CGNESolver.getId(), "CGNETestSolver" );

    CGNE CGNESolver2( "CGNETestSolver2" );
    BOOST_CHECK_EQUAL( CGNESolver2.getId(), "CGNETestSolver2" );

    CGNE CGNESolver3( CGNESolver2 );
    BOOST_CHECK_EQUAL( CGNESolver3.getId(), "CGNETestSolver2" );
    BOOST_CHECK( CGNESolver3.getPreconditioner() == 0 );

    CGNE CGNESolver4( "CGNESolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    CGNESolver4.setPreconditioner( preconditioner );

    CriterionPtr criterion( new IterationCount( 10 ) );
    CGNESolver4.setStoppingCriterion( criterion );

    CGNE CGNESolver5( CGNESolver4 );
    BOOST_CHECK_EQUAL( CGNESolver5.getId(), CGNESolver4.getId() );
    BOOST_CHECK_EQUAL( CGNESolver5.getPreconditioner()->getId(), CGNESolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
