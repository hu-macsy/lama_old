/**
 * @file CGSTest.cpp
 *
 * @license
 * Copyright (c) 2009-2016
 * Fraunhofer Institute for Algorithms and Scientific Computing SCAI
 * for Fraunhofer-Gesellschaft
 *
 * This file is part of the Library of Accelerated Math Applications (LAMA).
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
 * @endlicense
 *
 * @brief CGSTest.cpp
 * @author Jan Ecker
 * @date 09.03.2016
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
