/**
 * @file QMRTest.cpp
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
 * @brief QMRTest.cpp
 * @author lschubert
 * @date 07.08.2013
 */

#include <boost/test/unit_test.hpp>

#include <scai/solver/QMR.hpp>
#include <scai/solver/TrivialPreconditioner.hpp>
#include <scai/solver/criteria/IterationCount.hpp>
#include <scai/solver/logger/Timer.hpp>
#include <scai/solver/logger/CommonLogger.hpp>

#include <scai/solver/test/TestMacros.hpp>

using namespace scai::solver;

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE( QMRTest )

SCAI_LOG_DEF_LOGGER( logger, "Test.QMRTest" )

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE( ConstructorTest )
{
    LoggerPtr slogger( new CommonLogger( "<QMR>: ", LogLevel::noLogging, LoggerWriteBehaviour::toConsoleOnly ) );
    QMR QMRSolver( "QMRTestSolver", slogger );
    BOOST_CHECK_EQUAL( QMRSolver.getId(), "QMRTestSolver" );
    QMR QMRSolver2( "QMRTestSolver2" );
    BOOST_CHECK_EQUAL( QMRSolver2.getId(), "QMRTestSolver2" );
    QMR QMRSolver3( QMRSolver2 );
    BOOST_CHECK_EQUAL( QMRSolver3.getId(), "QMRTestSolver2" );
    BOOST_CHECK( QMRSolver3.getPreconditioner() == 0 );
    QMR QMRSolver4( "QMRSolver4" );
    SolverPtr preconditioner( new TrivialPreconditioner( "Trivial preconditioner" ) );
    QMRSolver4.setPreconditioner( preconditioner );
    CriterionPtr criterion( new IterationCount( 10 ) );
    QMRSolver4.setStoppingCriterion( criterion );
    QMR QMRSolver5( QMRSolver4 );
    BOOST_CHECK_EQUAL( QMRSolver5.getId(), QMRSolver4.getId() );
    BOOST_CHECK_EQUAL( QMRSolver5.getPreconditioner()->getId(), QMRSolver4.getPreconditioner()->getId() );
}

// ---------------------------------------------------------------------------------------------------------------

BOOST_AUTO_TEST_SUITE_END();
