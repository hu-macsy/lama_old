/**
 * @file RichardsonTest.cpp
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
 * @endlicense
 *
 * @brief Contains the implementation of the class RichardsonTest.
 * @author
 * @date 17.04.2015
 */
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
