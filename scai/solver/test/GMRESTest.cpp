/**
 * @file GMRESTest.cpp
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
 * @brief Contains the implementation of the class GMRESTest.
 * @author Malte Förster
 * @date 10.04.2012
 */

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
