/**
 * @file CGNRTest.cpp
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
 * @brief CGNRTest.cpp
 * @author Jan Ecker
 * @date 09.03.2016
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
